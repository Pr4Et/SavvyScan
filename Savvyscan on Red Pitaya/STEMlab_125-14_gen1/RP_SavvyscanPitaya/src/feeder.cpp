// feeder.cpp — C++20 feeder for Red Pitaya
// Designed by Shahar Seifer, Elbaum lab, Weizmann Institute of Science
// Based on original SavvyScan system
// Programming assisted by M365-Copilot

#include "feeder.h"
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>
#include <cmath>
#include <fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <sys/inotify.h>
#include <sys/select.h>
#include <cerrno>
#include <string>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <cstdarg>
#include <iostream>

#include <stdint.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <errno.h>

#include <sys/mman.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>


// Hardware headers
#include <rp.h>          // Red Pitaya C-API master header
#include "rp_gate.h"     // kernel driver's IOCTLs

using byte   = uint8_t; 
using uint64 = uint64_t;

#define MAX_PATERN_DATA_SIZE 0x00010000u // 64 KB
// --------- SoC base/size (UPDATE BASE_MEM to your Address Editor base) ---------
static constexpr uint32_t BASE_CTRL = 0x40600000u;  // sys[6]
static constexpr uint32_t SPAN_CTRL = 0x00001000u;
static constexpr uint32_t BASE_MEM  = 0x80000000u;  // <-- PUT YOUR BRAM base here
static constexpr uint32_t SPAN_MEM  = MAX_PATERN_DATA_SIZE;  
static constexpr uint32_t PAGESZ    = 4096u;
// --------- CTRL register offsets ----------
static constexpr uint32_t OFF_CTRL        = 0x00;
static constexpr uint32_t OFF_TICKDIV     = 0x04;
static constexpr uint32_t OFF_PWIDTH      = 0x08;
static constexpr uint32_t OFF_PLEN        = 0x0C;
static constexpr uint32_t OFF_MASK_COUNT  = 0x10;
static constexpr uint32_t OFF_MASK_STATE  = 0x14;
// CTRL bits
static constexpr uint32_t CTRL_ARM_BIT     = (1u << 0);
static constexpr uint32_t CTRL_STOP_BIT    = (1u << 1);
static constexpr uint32_t CTRL_SENDCAM_BIT = (1u << 2);
static constexpr uint32_t XLEN_ADDR  = 0x28;  // XLEN (pairs count) address in sys bus of the FPGA
static constexpr uint32_t OFF_TRIGGER_DELAY  = 0x2C;  

//Controlling the DMA that trasnfers DAC X,Y values from DDR to the FPGA:
static constexpr uint32_t DMA_BASE = 0x80400000u;   // DMA S_AXI_LITE base
static constexpr uint32_t DMA_SPAN = 0x00010000u;   // 64k is plenty for DMA regs

static int fd_mem =-1;
static volatile uint32_t* map_ctrl_mmap = (volatile uint32_t*) MAP_FAILED;
static volatile uint32_t* bram_mmap = (volatile uint32_t*) MAP_FAILED;
static volatile uint32_t* ddr_ptr_mmap = (volatile uint32_t*) MAP_FAILED;
static volatile uint32_t* dma_ctrl_mmap = (volatile uint32_t*) MAP_FAILED;
static volatile uint32_t* bram_mem32 = nullptr;
static volatile uint32_t* ctrl4mapctrl = nullptr;
static uint32_t dmm_size = 0;

static uint32_t pattern_number_of_active_entries=0; //number of 16-bit words we actuallt send, must be smaller than 32K

static int GUI_OutputAmpmV = 2000;
//static int GUI_BiasOutputP = -16; //F20
//static float output_amplifiers_gain_divider= 6.56;//F20
static int GUI_BiasOutputP = 0; //Krios
static float output_amplifiers_gain_divider= 4.4; //Krios
static float parking_voltage= 0.7;
static int ScanPattern =1;


static bool LowRes = false;
static int  LowResTime_uS = 13; //minimum duration between camera triggers, with zero delay, for which we treat the scan as 4d-stem
static bool use_ARINA = false;
static bool use_120kv = false;
static bool is_tilt_series = false;
static int  max_tilts = 1;
static std::string arinaFileName;
static int TiltIndex = 0;
static uint32_t CameraTrigger_delay  = 125;    //1 us @ 125 MHz

const std::string kMonitorPath = "/root/serialem_in.txt";

static int  recent_oversampling = 1;
static byte Gate_pattern[MAX_PATERN_DATA_SIZE]; 
static int  Gate_pattern_count = (int)(MAX_PATERN_DATA_SIZE);

static int  packet_id_last = 0;
static int EnergyKV=200;

static void* pvBufferREC;         // output image buffer
static void* pvBufferREC_buffer;  // per-line ADC buffer
static void* pvBufferREPxy;        // scan x,y table interveaved (DAC)
static void* pvBufferKey;         // (x,y) mapping table
static void* remote_data;

uint64 qwMemInBytesREC, qwMemInBytesREC_buffer;
uint64 qwMemInBytesREPxy;

static int samplerate, DACsamplerate_x,DACsamplerate_y;
static int NumDACx_samples, NumDACy_samples, NumDACxy_samples, ADCfetch_start, NumADC_samples;
static int NumSamplesPerCh, NumSamplesPerChREC, NumLocPerCh;

static bool has_output_reset= false;

const double PI = 3.141592653589793;
const int BytesPerSample = 2;
const int REPnumchannels = 2;
const int RECnumchannels = 1;

double MAT_fullXsize = 0;
double MAT_fullYsize = 0;
double MAT_netXsize = 0;
double MAT_netYsize = 0;
double MAT_netXstart = 0;
double MAT_netYstart = 0;
double MAT_samples_per_pixel = 0;
double MAT_exposure_us = 0;

//In the start bash script define your DECTRIS_SERVER_IP_ADDR="192.168.100.70" export DECTRIS_SERVER_IP_ADDR
const char* dectris_ip = getenv("DECTRIS_SERVER_IP_ADDR");

const char* home = getenv("HOME");

char msg[256];
static const uint8_t DIOP_BIT = 0; // DIO0_P


static size_t g_rec_bytes = 0, g_key_bytes = 0, g_repxy_bytes = 0;

static void free_scan_buffers() {
    if (pvBufferREC)  { delete[] reinterpret_cast<char*>(pvBufferREC);  pvBufferREC  = nullptr; g_rec_bytes  = 0; }
    if (pvBufferKey)  { delete[] reinterpret_cast<char*>(pvBufferKey);  pvBufferKey  = nullptr; g_key_bytes  = 0; }
    if (pvBufferREPxy) { delete[] reinterpret_cast<char*>(pvBufferREPxy); pvBufferREPxy = nullptr; g_repxy_bytes = 0; }
}


auto to_norm = [](float s) -> float {
    // Clamp hard to [-1.0, +1.0]
	return std::clamp(s, -1.0f, 1.0f);
};

inline int16_t to_dac(float v) {
    if (v > 1.0f) v = 1.0f;
    if (v < -1.0f) v = -1.0f;
    return (int16_t)std::lrint(v * 8191.0f);
}




static inline void LogFile(const char* fmt, ...) {
    FILE* fp = std::fopen("/tmp/savvyscan-feeder.log", "a");
    if (!fp) return;
    va_list ap; va_start(ap, fmt);
    std::vfprintf(fp, fmt, ap);
    std::fprintf(fp, "\n");
    std::fflush(fp);
    va_end(ap);
    std::fclose(fp);
}

// ~ tilde expansion helpers / inotify monitor (unchanged) --------------------
static std::string expand_tilde(const std::string& p) {
    if (!p.empty() && p[0] == '~') {
        const char* h = std::getenv("HOME"); if (!h) h = "/";
        return (p.size() == 1) ? std::string(h) : std::string(h) + p.substr(1);
    }
    return p;
}

static int g_inotify_fd = -1;
static int g_watch_fd   = -1;
static std::string g_mon_path;

bool monitor_init() {
    g_mon_path = expand_tilde(kMonitorPath);
    g_inotify_fd = inotify_init1(IN_CLOEXEC);
    if (g_inotify_fd < 0) { perror("inotify_init1"); return false; }
    uint32_t mask = IN_CLOSE_WRITE | IN_MODIFY | IN_MOVED_TO;
    g_watch_fd = inotify_add_watch(g_inotify_fd, g_mon_path.c_str(), mask);
    if (g_watch_fd < 0) { perror("inotify_add_watch"); close(g_inotify_fd); g_inotify_fd = -1; return false; }
    std::printf("[monitor] watching %s\n", g_mon_path.c_str());
    return true;
}

void monitor_shutdown() {
    if (g_watch_fd >= 0) { inotify_rm_watch(g_inotify_fd, g_watch_fd); g_watch_fd = -1; }
    if (g_inotify_fd >= 0) { close(g_inotify_fd); g_inotify_fd = -1; }
}

static std::string trim(std::string s) {
    // Lambda predicate: std::isspace expects unsigned char (or EOF)
    auto isspace_uc = [](unsigned char ch){ return std::isspace(ch); };
    // left trim: find first non-space
    auto it = s.begin();
    while (it != s.end() && isspace_uc(static_cast<unsigned char>(*it))) {
        ++it;
    }
    s.erase(s.begin(), it);
    // Right trim loop
    auto rit = s.rbegin();
    while (rit != s.rend() && isspace_uc(static_cast<unsigned char>(*rit))) {
        ++rit;
    }
    s.erase(rit.base(), s.end());
    return s;
}

static bool starts_with(const std::string& s, const char* prefix) {
    return s.size() >= std::strlen(prefix) &&
           std::equal(prefix, prefix + std::strlen(prefix), s.begin());
}



static void apply_command_line(const std::string& line_in)
{

    std::string line = trim(line_in);
    if (line.empty()) return;
	std::fprintf(stderr, "%s\n", line_in.c_str());
    if (starts_with(line, "ARINA:")) {
        std::string fname = trim(line.substr(6));
        use_ARINA = !fname.empty();
        arinaFileName = fname;
        std::printf("[monitor] ARINA -> %s, use_ARINA=%d\n", arinaFileName.c_str(), use_ARINA);
        return;
    }
    if (starts_with(line, "StartTiltSeries:")) {
        std::string num = trim(line.substr(std::strlen("StartTiltSeries:")));
        int n = 0; try { n = std::stoi(num); } catch (...) { n = 0; }
        if (n <= 0) n = 1;
        is_tilt_series = true; TiltIndex = 0; max_tilts = n;
        std::printf("[monitor] StartTiltSeries -> max_tilts=%d, is_tilt_series=1\n", max_tilts);
        return;
    }
    if (starts_with(line, "SetThTimeus:")) {
        std::string num = trim(line.substr(std::strlen("SetThTime:")));
        int n = 0; try { n = std::stoi(num); } catch (...) { n = 0; }
        if (n <= 0) n = 1;
        LowResTime_uS = n;
        std::printf("[monitor] SetThTime us -> %d \n", LowResTime_uS);
        return;
    }
	if (starts_with(line, "SetTriggerDelayus:")) {
        std::string num = trim(line.substr(std::strlen("SetTriggerDelayus:")));
        int n = 0; try { n = std::stoi(num); } catch (...) { n = 0; }
        if (n <= 0) n = 1;
        CameraTrigger_delay = n*125;
        std::printf("[monitor] SetTriggerDelayus -> %d \n", CameraTrigger_delay/125);
        return;
    }

    if (line == "StopTiltSeries") {
        is_tilt_series = false; TiltIndex=0; use_ARINA=false;
        std::printf("[monitor] StopTiltSeries -> is_tilt_series=0\n");
        return;
    }
    if (starts_with(line, "EnergyKV:")) {
        std::string num = trim(line.substr(std::strlen("EnergyKV:")));
        int n = 0; try { n = std::stoi(num); } catch (...) { n = 0; } 
        if (n>0) EnergyKV=n;
        std::printf("[monitor] EnergyKV -> %d \n", EnergyKV);
        return;
    }
    if (starts_with(line, "ScanPattern:")) {
        std::string num = trim(line.substr(std::strlen("ScanPattern:")));
        int n = 0; try { n = std::stoi(num); } catch (...) { n = 0; }
        if (n <= 0) n = 1;
        ScanPattern = n;
        std::printf("[monitor] ScanPattern -> %d \n", ScanPattern);
        return;
    }
	
}

static void read_and_apply_all_lines() {
    std::ifstream f(g_mon_path);
    if (!f) { perror("[monitor] open"); return; }
    std::string line;
    while (std::getline(f, line)) apply_command_line(line);
}

void monitor_pump_once() {
    if (g_inotify_fd < 0) return;
    fd_set rfds; FD_ZERO(&rfds); FD_SET(g_inotify_fd, &rfds);
    timeval tv{0,0};
    int rv = select(g_inotify_fd + 1, &rfds, nullptr, nullptr, &tv);
    if (rv <= 0) return;
    char buf[8192];
    ssize_t n = ::read(g_inotify_fd, buf, sizeof(buf));
    if (n < 0) { if (errno != EAGAIN && errno != EINTR) perror("inotify read"); return; }
    for (char* p = buf; p < buf + n; ) {
        auto* ev = reinterpret_cast<inotify_event*>(p);
        if (ev->wd == g_watch_fd && (ev->mask & (IN_CLOSE_WRITE | IN_MOVED_TO | IN_MODIFY))) {
            read_and_apply_all_lines();
        }
        p += sizeof(inotify_event) + ev->len;
    }
}

// Pulse DIO0_P for 'pulse_us' microseconds
static bool pulse_dio0p_us(unsigned pulse_us)
{
    if (rp_DpinSetDirection(RP_DIO0_P, RP_OUT) != RP_OK) return false;
    if (rp_DpinSetState(RP_DIO0_P, RP_HIGH)     != RP_OK) return false;
    usleep(pulse_us);
    if (rp_DpinSetState(RP_DIO0_P, RP_LOW)      != RP_OK) return false;
    return true;
}

// Quantize samplerate to ADC decimation factor (1..65536)
static uint32_t pick_decimation(double samplerate_hz) {
    double target = samplerate_hz;
    double best_err = 1e300;
    uint32_t best_dec = 1;
    for (uint32_t dec = 1; dec <= 65536; ++dec) {
        double rate = 125e6 / (double)dec;
        double err  = std::abs(rate - target);
        if (err < best_err) { best_err = err; best_dec = dec; }
    }
    return best_dec;
}


// Align 'x' upward to the next multiple of 4096 bytes.
static inline uint32_t ALIGN_TO_4096(uint32_t x)
{
    const uint32_t PAGE = 4096u;
    return (x + PAGE - 1u) & ~(PAGE - 1u);
}

// Align down to page for mmap
static inline off_t page_align_down(uint32_t addr) {
    long pagesz = sysconf(_SC_PAGESIZE);
    return (off_t)(addr & ~(pagesz - 1));
}

int copy_dmm_to_adc16(uint32_t adc_buffer_start, uint32_t dmm_start,
                      uint32_t words_to_copy,
                      int16_t *dst)
{
    // Calculate offset inside the already-mapped DDR region
    if (adc_buffer_start < dmm_start) {
        fprintf(stderr, "ADC buffer start below DDR base\n");
        return -1;
    }

    uint32_t byte_offset = adc_buffer_start - dmm_start;

    // Invalidate cache BEFORE reading DMA-written DDR
    msync((void *)((uint8_t *)ddr_ptr_mmap + byte_offset),
          words_to_copy * 4,
          MS_INVALIDATE);

    const volatile uint32_t *src32 =
        ddr_ptr_mmap + (byte_offset / 4);

    // Unpack: 2 x 16-bit signed samples per 32-bit word
    for (uint32_t i = 0; i < words_to_copy; ++i) {
        uint32_t w = src32[i];

        dst[2*i + 0] = (int16_t)( w        & 0xFFFF);
        dst[2*i + 1] = (int16_t)((w >> 16) & 0xFFFF);
    }

    return 0;
}



static inline uint32_t reg_rd(volatile uint32_t* regs, uint32_t off) {
    return regs[off/4];
}
static inline void reg_wr(volatile uint32_t* regs, uint32_t off, uint32_t val) {
    regs[off/4] = val;
}


// --------- mmap helper ----------
static uint32_t* map_region32(int fd,
                              uint32_t base,
                              uint32_t span,
                              volatile uint32_t** out32)
{
    long pagesz = sysconf(_SC_PAGESIZE);
    off_t page_base = base & ~(off_t)(pagesz - 1);
    off_t page_off  = base - page_base;
    // Enforce 32-bit alignment
    if ((base & 0x3) != 0) {
        std::fprintf(stderr,
                     "map_region32: base 0x%08X not 32-bit aligned\n",
                     base);
        return nullptr;
    }
    // mmap return is untyped; cast once
    uint32_t* map = static_cast<uint32_t*>(
        mmap(nullptr,
             span + page_off,
             PROT_READ | PROT_WRITE,
             MAP_SHARED,
             fd,
             page_base)
    );
    if (map == MAP_FAILED) {
        std::fprintf(stderr,
                     "mmap(0x%08X, 0x%X) failed: %s\n",
                     base, span, std::strerror(errno));
        return nullptr;
    }
    // Offset in BYTES, then recast as 32-bit
    auto* regs32 = reinterpret_cast<volatile uint32_t*>(
        reinterpret_cast<uint8_t*>(map) + page_off
    );
    *out32 = regs32;
    return map;
}



static bool write_table_to_ddr(volatile uint32_t* ddr,
                               uint32_t offset_bytes,
                               const int16_t* table,
                               size_t count)
{
	volatile uint8_t* dst = (volatile uint8_t*)ddr + offset_bytes;
    memcpy((void*)dst, table, count * sizeof(int16_t));
    return true;
}



static const double FMIN_SET = 0.1; // lower bound for configured FREQ stability

uint32_t choose_scale_Ny(uint32_t Ny, double F_frame_phys) {
    if (Ny == 0 || F_frame_phys <= 0.0) return 1;  // safety
    // k to satisfy the frequency floor
    double k_freq = std::ceil((FMIN_SET * 16384.0) / (F_frame_phys * (double)Ny));
    // k to reach/cap the 16384 samples
    double k_cap  = std::ceil(16384.0 / (double)Ny);
    // Minimum integer that satisfies "either" condition:
    uint32_t scale_Ny = (uint32_t)std::max(1.0, std::min(k_freq, k_cap));
    return scale_Ny;
}


// Helper: set AWG frequencies and read back actual implemented values
struct AwgFreqs {
    double F1_req; // requested CH1 frequency
    double F2_req; // requested CH2 frequency
    double F1_act; // actual implemented CH1 frequency
    double F2_act; // actual implemented CH2 frequency
};



static inline bool dio0p_is_high() {
    rp_pinState_t s = RP_LOW;
    rp_DpinGetState(RP_DIO0_P, &s);
    return s == RP_HIGH;
}

static void wait_dio0p_low() {
    // Ensure we see a low period before arming (avoid arming inside a HIGH pulse)
    for (int i = 0; i < 10000; ++i) { // ~10 ms worst-case guard
        if (!dio0p_is_high()) return;
        usleep(1);
    }
}

// -------------------------------------------------------------
// Start AXI DMA MM2S to stream XY data from DDR to FPGA
// base_addr: physical DDR address of first XY pair
// num_pairs: number of XY pairs (each 4 bytes: {X,Y})
// -------------------------------------------------------------
static void start_dma_mm2s(uint32_t base_addr, uint32_t num_pairs)
{
	

    //static constexpr uint32_t DMA_BASE = 0x80400000u;   // your BD address
    //static constexpr uint32_t DMA_SPAN = 0x00010000u;   // enough for dma_ctrl_mmap regs

	dma_ctrl_mmap = (volatile uint32_t*)
		mmap(NULL, DMA_SPAN,
			 PROT_READ | PROT_WRITE,
			 MAP_SHARED,
			 fd_mem,
			 DMA_BASE);
	if (dma_ctrl_mmap == MAP_FAILED) {
		perror("mmap DMA failed");
		return;
	}
	
    const uint32_t MM2S_DMACR   = 0x00/4;
	const uint32_t MM2S_DMASR   = 0x04/4;
    const uint32_t MM2S_SA      = 0x18/4;
    const uint32_t MM2S_LENGTH  = 0x28/4;

	const uint32_t DMACR_RUN = 1 << 0;
	const uint32_t DMACR_RESET = 1 << 2;
	const uint32_t DMACR_IOC = 1 << 12;
	const uint32_t DMACR_ERR = 1 << 14;

    // Reset DMA channel
    dma_ctrl_mmap[MM2S_DMACR] = DMACR_RESET;   // Reset by setting RESET bit
    usleep(5);
	while (dma_ctrl_mmap[MM2S_DMACR] & DMACR_RESET) {}   // wait reset complete
	//printf("After Reset: MM2S_DMASR=0x%08X, DMACR=0x%08X, DMA MM2S_LENGTH=%u, DMA MM2S_SA=0x%08X\n",dma_ctrl_mmap[MM2S_DMASR],dma_ctrl_mmap[MM2S_DMACR],dma_ctrl_mmap[MM2S_LENGTH],dma_ctrl_mmap[MM2S_SA]);
	// Set source
	dma_ctrl_mmap[MM2S_SA] = base_addr;  // must be aligned
	// Start engine
	dma_ctrl_mmap[MM2S_DMACR] = DMACR_RUN | DMACR_IOC | DMACR_ERR;//  RS=1, IOC_IrqEn=1, Err_IrqEn=1
	// Start transfer
	dma_ctrl_mmap[MM2S_LENGTH] = num_pairs * 4;
	
	const uint32_t TIMEOUT_US = 30000;  // 30 ms
	uint32_t waited = 0;
	while (true) {
		uint32_t sr = dma_ctrl_mmap[MM2S_DMASR];
		if ((sr & (1<<12)) || (sr & (1<<1))) break;     // IOC – done or idle
		if (sr & (1<<14)) {          // Error
			printf("DMA ERROR: DMASR=0x%08X\n", sr);
			break;
		}
		if (waited >= TIMEOUT_US) {
			//printf("DMA uploading finished by timeout\n");
			break;
		}
		usleep(1);
		waited++;
	}
	printf("After transfer: MM2S_DMASR=0x%08X, DMACR=0x%08X, DMA MM2S_LENGTH=%u, DMA MM2S_SA=0x%08X, num_pairs=%u\n",dma_ctrl_mmap[MM2S_DMASR],dma_ctrl_mmap[MM2S_DMACR],dma_ctrl_mmap[MM2S_LENGTH],dma_ctrl_mmap[MM2S_SA], num_pairs);
	
}




static inline void reg_wr32(volatile uint32_t* regs, uint32_t off, uint32_t val) {
    regs[off/4] = val;
}
static inline uint32_t reg_rd32(volatile uint32_t* regs, uint32_t off) {
    return regs[off/4];
}

// --------- Pattern writer: 16-bit entries → 32-bit words (RMW-safe) ---------
// Writes 'entries' 16-bit values into 64KB BRAM window.
// BRAM word i contains entries [2*i] (low half) and [2*i+1] (high half).
static bool write_pattern_16(volatile uint32_t* tmem32,
                             const int16_t* pat16,
                             size_t entries)
{
    if (entries > 32768) { // 64KB / 2B
        std::fprintf(stderr, "Pattern too large: %zu (max 32768 16-bit entries)\n", entries);
        return false;
    }
    // Write in 32-bit words using read-modify-write, so we never clobber the other half.
    size_t full_words = entries / 2;
    size_t tail       = entries & 1u;

    // Write all full 32-bit words (two 16-bit entries per word)
    for (size_t w = 0; w < full_words; ++w) {
        uint16_t lo = static_cast<uint16_t>(pat16[2*w + 0]);
        uint16_t hi = static_cast<uint16_t>(pat16[2*w + 1]);
        uint32_t packed = (static_cast<uint32_t>(hi) << 16) | lo;
        tmem32[w] = packed;
    }

    // If an odd entry remains, RMW the last word's upper/lower half
    if (tail) {
        size_t w = full_words;
        uint32_t cur = tmem32[w];
        uint16_t last = static_cast<uint16_t>(pat16[2*w]); // remaining entry
        // Put it in LOW half (by convention we fill low first). If you want it to be high half,
        // change the other half accordingly.
        uint32_t newv = (cur & 0xFFFF0000u) | last;
        tmem32[w] = newv;
    }
    return true;
}


int set_rp_gate_and_acqstart(uint32_t div,bool send_camera_trigger,uint32_t dmm_start,uint32_t XY_offset)
{
   // Cast to 16-bit entries
    const int16_t* Gate_pattern_words = reinterpret_cast<const int16_t*>(Gate_pattern);
    const int16_t* AWGtablexy = static_cast<const int16_t*>(pvBufferREPxy);
	
    // Map CTRL and MEM
    map_ctrl_mmap = map_region32(fd_mem, BASE_CTRL, SPAN_CTRL, &ctrl4mapctrl);
    if (!map_ctrl_mmap) { return 2; }
	
    bram_mmap = map_region32(fd_mem, BASE_MEM,  SPAN_MEM,  &bram_mem32);  //returned address: bram_mem32
    if (!bram_mmap)  { return 3; }

    // --------- Configure control registers ----------
    const uint32_t TICKDIV = div;   // dividing 125MHz clock to ADC clock
    const uint32_t PWIDTH  = 375;    //3 us @ 125 MHz
    const uint16_t MASK_ST = 0x8000;   // MSB = state
    const uint16_t MASK_CT = 0x7FFF;   // lower 15 bits = count

    // Clear CTRL
    reg_wr32(ctrl4mapctrl, OFF_CTRL, 0);

	// Number of XY PAIRS (each pair = {X,Y})
	uint32_t NumXYpairs = NumDACxy_samples / 2;
	// Tell FPGA the DDR memory size for DAC 
	reg_wr32(ctrl4mapctrl, XLEN_ADDR,  NumXYpairs);                                                                                                                                                                                                                                                                                                                                   
	//No beam parking starting here
    reg_wr32(ctrl4mapctrl, 0x1C, 0);      // PARK_EN = 0, , no parking of beam  


    // Program masks/timing
    reg_wr32(ctrl4mapctrl, OFF_MASK_STATE, MASK_ST);
    reg_wr32(ctrl4mapctrl, OFF_MASK_COUNT, MASK_CT);
    reg_wr32(ctrl4mapctrl, OFF_TICKDIV,    TICKDIV);
    reg_wr32(ctrl4mapctrl, OFF_PWIDTH,     PWIDTH);
    reg_wr32(ctrl4mapctrl, OFF_PLEN, pattern_number_of_active_entries);
	reg_wr32(ctrl4mapctrl, OFF_TRIGGER_DELAY,     CameraTrigger_delay);

    // --------- Write the pattern into BRAM ----------
    // Pack two 16-bit entries per 32-bit word, RMW-safe on odd tail
    if (!write_pattern_16(bram_mem32, Gate_pattern_words, pattern_number_of_active_entries)) {
        std::fprintf(stderr, "Pattern write failed\n");
        return 4;
    }

	start_dma_mm2s(dmm_start + XY_offset, NumXYpairs);
	
	
	//request acquisition triggered by a pulse, but expect it will take affect only after the first arm pulse that starts AWG alone.
	rp_AcqStart();  //movedbefore arm pulse

    // --------- ARM sequence ----------
    uint32_t ctrl_val = send_camera_trigger ? CTRL_SENDCAM_BIT : 0u;

    // ARM rising edge
    reg_wr32(ctrl4mapctrl, OFF_CTRL, ctrl_val | CTRL_ARM_BIT);
    usleep(3); //3us — allow your ARM detector to see the rising edge
    reg_wr32(ctrl4mapctrl, OFF_CTRL, ctrl_val & ~CTRL_ARM_BIT);

    std::fprintf(stderr, "Pattern loaded: %u entries (16-bit). LowRes=%d\n", pattern_number_of_active_entries,LowRes);

 
	return 0;
}

int stop_rp_gate() {


    // ---------------------------------------------------
    // 1. Enable parking mode
    // ---------------------------------------------------
    reg_wr32(ctrl4mapctrl, 0x1C, 1);                        // PARK_EN = 1
    reg_wr32(ctrl4mapctrl, 0x20, to_dac(parking_voltage));  // PARK_X
    reg_wr32(ctrl4mapctrl, 0x24, to_dac(parking_voltage));  // PARK_Y

    // ---------------------------------------------------
    // 2. Optionally stop the FSM
    // ---------------------------------------------------
    reg_wr32(ctrl4mapctrl, OFF_CTRL, CTRL_STOP_BIT);


    return 0;
}


bool run_scan_with_gate_and_adc()
{
    const int16_t* AWGtablexy = static_cast<const int16_t*>(pvBufferREPxy); //previously float float
    int16_t*     ADCtable  = static_cast<int16_t*>(pvBufferREC);

    // ---- 1) Initialize Red Pitaya C API ----
    if (rp_Init() != RP_OK) return false;

	
    // ---- 2) Configure acquisition (ADC) ----
    uint32_t dec    = pick_decimation(samplerate);
    double   Fs_ADC = 125e6 / (double)dec;   // exact in FPGA decimator
	
	const double F_frame_phys = Fs_ADC  / (double)NumDACy_samples;       // physcial frames per second

    rp_AcqReset();
	uint32_t number_of_32bit_words = qwMemInBytesREC / 4;
	uint32_t qwMemInPointsREC = number_of_32bit_words * 2;	

	//-----------------------------------------------------
	// Map DDR region (32 MB deep memory acquisition region)
	//-----------------------------------------------------
	uint32_t dmm_start = 0;
	rp_AcqAxiGetMemoryRegion(&dmm_start, &dmm_size);

	ddr_ptr_mmap = (volatile uint32_t*)
		mmap(NULL, dmm_size,
			 PROT_READ | PROT_WRITE,
			 MAP_SHARED,
			 fd_mem,
			 dmm_start);
	if (ddr_ptr_mmap == MAP_FAILED) {
		perror("mmap DDR failed");
		return false;
	}
	//-----------------------------------------------------
	// 1. Write XY DAC table (interleaved) into DDR
	//-----------------------------------------------------
	uint32_t XY_offset = 0;   // table starts at beginning of region
	uint32_t XY_bytes  = NumDACxy_samples * sizeof(int16_t);
	write_table_to_ddr(ddr_ptr_mmap, XY_offset,
					   (const int16_t*)AWGtablexy,
					   NumDACxy_samples);
	//-----------------------------------------------------
	// 2. Place ADC capture buffer AFTER the XY table
	//-----------------------------------------------------
	uint32_t ADC_offset = ALIGN_TO_4096(XY_bytes);
	uint32_t adc_buffer_start = dmm_start + ADC_offset;
	// Align ADC capture to a 4096B boundary (required by DMA)
	uint32_t adc_words = qwMemInBytesREC / 4;            // number of 32-bit DMA words
	uint32_t adc_bytes = adc_words * 4;
	uint32_t adc_bytes_aligned = ALIGN_TO_4096(adc_bytes);
	uint32_t adc_samples       = adc_bytes_aligned / 2;  // 2 bytes per 16-bit sample
	// Program deep memory acquisition region
	rp_AcqAxiSetBufferSamples(RP_CH_1,
							  adc_buffer_start,
							  adc_samples);
	
		
	// Usual DMM configuration (decimation, trigger delay, etc.)
	rp_AcqAxiSetDecimationFactor(dec);
	rp_AcqAxiSetTriggerDelay(RP_CH_1, 0);
	// Set global trigger to external DIO0_P rising edge
	rp_AcqSetTriggerSrc(RP_TRIG_SRC_EXT_PE);
	// Enable and run acquisition
	rp_AcqAxiEnable(RP_CH_1, true);
	
	fprintf(stderr,
			"dec=%u  Fs_ADC=%.3f kS/s  N_positions=%d, ADCfetch_start=%d\n",
			dec, Fs_ADC/1e3, NumDACx_samples, ADCfetch_start);

    bool send_camera_trigger =!LowRes; //(!LowRes && use_ARINA);
	
	usleep(250000); //wait 0.25 sec, since otherwise the initial triggers are quiet
	
	int answer=set_rp_gate_and_acqstart(dec, send_camera_trigger,dmm_start,XY_offset);
	if (answer!=0)
	{
		fprintf(stderr,"Error in setting FPGA. Err num= %u",answer);
		return false;
	}
	
	//RUN mass acquisition
	usleep((uint32_t)(1000000.0/F_frame_phys+1000000.0/(double)samplerate)+1); //wait in microseconds
	// ------------------- WAIT FOR FILL -----------------------------


	// ------------------- COPY DATA OUT OF DDR ----------------------
	if (copy_dmm_to_adc16(adc_buffer_start, dmm_start, number_of_32bit_words, ADCtable) != 0) {
		fprintf(stderr, "Error copying DMM data\n");
		return false;
	}
	usleep(5000);
	stop_rp_gate();
	//-----------------------------------------------------


	usleep(5000);    //starngely needs a delay, but do not exceed becuase the ADC table could be overwritten 


    rp_Release();
	usleep(5000);
    return true;
}

static void msleep(unsigned ms) {
    timespec ts{ static_cast<time_t>(ms / 1000), static_cast<long>((ms % 1000) * 1000000) };
    nanosleep(&ts, nullptr);
}

// (Other functions below unchanged: send_pattern, run_Arina_acquire, prepare_AWG, correctedLOCxy, acquire_image, activate_scan, etc.)

bool send_pattern(const uint8_t* buf, size_t bytes, int dio_p_bit = 2)
{
    int fd = ::open("/dev/rp_gate", O_WRONLY);
    if (fd < 0) { perror("open"); return false; }
    rp_gate_cfg cfg{};
    cfg.mask_count          = 0x7FFF;
    cfg.mask_state          = 0x8000;
    cfg.dio_p_bit           = (uint8_t)dio_p_bit;
    cfg.start_bit           = 0;
    cfg.camera_bit          = 3;
    cfg.adcstart_bit        = 1;
    cfg.send_camera_trigger = 1;
    if (ioctl(fd, RP_GATE_CONFIG, &cfg) < 0) { perror("ioctl CONFIG"); ::close(fd); return false; }
    ssize_t wr = ::write(fd, buf, bytes);
    if (wr != (ssize_t)bytes) { perror("write"); ::close(fd); return false; }
    if (ioctl(fd, RP_GATE_ARM) < 0) { perror("ioctl ARM"); ::close(fd); return false; }
    ::close(fd);
    return true;
}



void run_Arina_acquire(double integration_time_us, int number_of_triggers)
{
    double t_sec = (integration_time_us - 3.0-(double)CameraTrigger_delay/125.0) * 1e-6;
    if (t_sec < 0.0) t_sec = 0.0;
 
    const char* python_exec = "/usr/bin/python3";
	if (!dectris_ip) {
		dectris_ip = "192.168.100.70"; // fallback if not set
	}


    char cmd[512];

    std::snprintf(
        cmd, sizeof(cmd),
        "%s -u -m acquireimages -i %s -x -t %.6f -n %d -e %d -o %s_s%d",
        python_exec, dectris_ip, t_sec, number_of_triggers, EnergyKV,
        arinaFileName.c_str(), TiltIndex
    );

    // Prepare ~/arm.sh content
    std::string scriptPath = std::string(home ? home : ".") + "/arm.sh";

    std::string script =
        "#!/bin/sh\n"
        "export PYTHONPATH=$PYTHONPATH:/opt/velan/python:/opt/velan/Lib\n"
        + std::string(cmd) + "\n";

    // Write script
    {
        FILE* f = fopen(scriptPath.c_str(), "w");
        if (!f) {
            fprintf(stderr, "Failed to write %s\n", scriptPath.c_str());
            return;
        }
        fwrite(script.data(), 1, script.size(), f);
        fclose(f);
        chmod(scriptPath.c_str(), 0755);
    }

    // Fork — child will exec the script, printing to same terminal
    pid_t pid = fork();
    if (pid < 0) {
        perror("fork");
        return;
    }
    if (pid == 0) {
        // CHILD: no redirection, inherit stdout/stderr
        execl(scriptPath.c_str(), scriptPath.c_str(), (char*)nullptr);

        // Only reached if exec fails
        perror("exec");
        _exit(127);
    }

    // PARENT: just continues
    printf("Launched acquisition script as PID %d\n", pid);
}



float checkLOCxy(float result)
{
	if (result > 32766.0) result = 32766.0;
	if (result < -32767.0) result = -32767.0;
	return result;

}

void checkworking() 
{
	LogFile("working");
}

/*
**************************************************************************
Prepare AWG (Arbitraty wave generator) data for scanx/y 
**************************************************************************
*/
void prepare_AWG(int binning, int scXsize, int scYsize, double AspectRatio, double rotation, int oversampling, double pixeltime_us,int ScanMode,double flyback_us, int width, int height,double sampling_pixeltime_us, int fsize) 
{
	
	//set_cards_samplerate_andmore(double sampling_pixeltime_us, int fsize, int& oversampling)
	oversampling=1;
	NumSamplesPerCh = fsize * oversampling;
	NumLocPerCh = fsize;
	
	int flypix = (int)(flyback_us / pixeltime_us);
	int samples_per_line = scXsize * oversampling + flypix * oversampling;
	int locs_per_line    = scXsize + flypix;

	
	
	DACsamplerate_x = (int)(1000000.0 / sampling_pixeltime_us); //For scanning pattern we only ask one point per pixel
	samplerate = DACsamplerate_x * oversampling;//ADC clock
    qwMemInBytesREC = NumSamplesPerCh * BytesPerSample * RECnumchannels *oversampling;
	qwMemInBytesREPxy = scYsize * locs_per_line * sizeof(int16_t)*2; 
	
	
	//const size_t need_rep  = qwMemInBytesREP;                
	const size_t need_rec  = qwMemInBytesREC;                      // int16 samples buffer
	const size_t need_key  = size_t(NumSamplesPerCh) * 2 * sizeof(int16_t);  // locx,locy as int16_t
	const size_t need_repxy = qwMemInBytesREPxy;                     // size of int16 full table, interleaved x,y values

	if (need_rec != g_rec_bytes) {
		if (pvBufferREC) delete[] reinterpret_cast<char*>(pvBufferREC);
		pvBufferREC  = static_cast<void*>(new char[need_rec]);
		g_rec_bytes  = need_rec;
	}
	if (need_key != g_key_bytes) {
		if (pvBufferKey) delete[] reinterpret_cast<char*>(pvBufferKey);
		pvBufferKey  = static_cast<void*>(new char[need_key]);
		g_key_bytes  = need_key;
	}
	if (need_repxy != g_repxy_bytes) {
		if (pvBufferREPxy) delete[] reinterpret_cast<char*>(pvBufferREPxy);
		pvBufferREPxy = static_cast<void*>(new char[need_repxy]);
		g_repxy_bytes = need_repxy;
	}

	
	NumDACx_samples=scYsize * locs_per_line; //previously locs_per_line;
	NumDACy_samples=scYsize * locs_per_line; //previously scYsize;
	NumDACxy_samples=scYsize * locs_per_line *2; //previously scYsize;
	NumADC_samples=scYsize*width; //not used, just informative 
	
	
	//fprintf(stderr, "Requested pixeltime_us=%.6f s, sampling_pixeltime_us=%.3f, fsize=%d \n",  pixeltime_us, sampling_pixeltime_us, fsize);

	LowRes=((int)((double)1000000/(double)samplerate -((double) CameraTrigger_delay/125.0)) < LowResTime_uS);
	
	double standard_amplitude=GUI_OutputAmpmV;
	int xscan_millivolts_amp=int(standard_amplitude*(scXsize*binning)/8192.0);
	int yscan_millivolts_amp=int(standard_amplitude* AspectRatio * (scYsize*binning)/8192.0);
	double orientation = 0;
	int xscan_millivolts_amp_tag = xscan_millivolts_amp;
	int yscan_millivolts_amp_tag = yscan_millivolts_amp;
	double scanx,scany;//would be scan in millivolt if not for rotation
	float DACx_tag=0;
	float DACy_tag=0;
	double prevLOCx=0;
	double prevLOCy=0;
	int16_t locx;
	int16_t locy;
	float cor_locx=0.0; 
	float cor_locy=0.0; 
	int marginx,marginy;
	int leftover=0;
	int counter=0;
	int counterLoc = 0;
	int rnbs,phnmbs;
	int scanr;
	double scanph,scanphcor;
	double xdrift;
	int epoch = 0;
	double Tcx = 0;
	double Tcy = 0;
	double sintheta=0;
	double costheta=1.0;
	//double delay=0.000220; //FEI deflection coils L/R, in seconds 
	marginx=(int)(0.5*(scXsize-width));
	marginy=(int)(0.5*(scYsize-height));
	
	MAT_fullXsize = scXsize;
	MAT_fullYsize = scYsize;
	MAT_netXsize = width;
	MAT_netYsize = height;
	MAT_netXstart = marginx;
	MAT_netYstart = marginy;
	MAT_exposure_us = pixeltime_us;


	
	//The numbers are stored in short format (signed 16bits)
	int16_t * AWGtablexy=(int16_t *)pvBufferREPxy; //previously float float
	int16_t * Keytable=(int16_t *)pvBufferKey;
	
	//here is scan mode 1. Doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left
	ADCfetch_start=0;
	float BiasOutput_mv = 0; 
	if (ScanPattern == 1) BiasOutput_mv = ((float)standard_amplitude * 0.414*(1.0+ (float)GUI_BiasOutputP/100.0) );

	switch (ScanPattern)
	{
	case 1: //### standard raster scan: Doing: scanx(line) sweeps from left (2^15-1) to right (-2^15) and then flyback over a period flyback_us back to most left
		for (int scanyvar = scYsize-1; scanyvar >=0; scanyvar--) //order reversed 8May2022
		{
			//scan from left to right
			scany = (float)yscan_millivolts_amp * ((float)scanyvar - (float)scYsize / 2.0) / ((float)scYsize / 2.0);
			DACy_tag = 0.001*scany/output_amplifiers_gain_divider; //(int)(scany * 32768.0 / (double)yscan_millivolts_amp_tag);
			for (int scanxvar = 0; scanxvar < scXsize * oversampling; scanxvar++)
			{
				scanx = (float)xscan_millivolts_amp * ((float)scanxvar - 0.5*(float)scXsize * (float)oversampling) / ((float)(scXsize * oversampling) / 2.0);
				DACx_tag= 0.001*(-scanx+BiasOutput_mv)/output_amplifiers_gain_divider;
				if (counter % oversampling == 0)
				{
					{
						*AWGtablexy++ = to_dac(DACx_tag);//store scanx_tag value
						*AWGtablexy++ = to_dac(DACy_tag);//store scany_tag value
					}
					counterLoc++;
				}
				locx = (int16_t)(scanxvar / oversampling - 2 * marginx);
				locy = (int16_t)(scanyvar - marginy);
			
				if (ADCfetch_start==0 && ((int)locx >= 0 && (int)locy >= 0 && (int)locx < width && (int)locy < height))
				{
					ADCfetch_start=counterLoc-1;
				}

				
				//cor_locx=correctedLOCxy(locx, delay, samplerate, counter, & prevLOCx);  
				//cor_locy=correctedLOCxy(locy, delay, samplerate, counter, & prevLOCy);
				//cor_locx = checkLOCxy(locx);
				//cor_locy = checkLOCxy(locy);
				counter++;
				*Keytable++ = locx; //location in the final image (+margin) for the data to fit into
				*Keytable++ = locy;
			}
			//flyback
			//int flypix =(int)(flyback_us / pixeltime_us);
			for (int scanxvar = flypix * oversampling - 1; scanxvar >= 0; scanxvar--)
			{
				scanx = xscan_millivolts_amp* (float(width + 2 * marginx) * (float)scanxvar / (float)flypix + ((-0.5 * (float)width - (float)marginx) * (float)oversampling)) / ((float)(scXsize * oversampling) / 2.0);
				DACx_tag=0.001*(-scanx+BiasOutput_mv)/output_amplifiers_gain_divider;
				if (counter % oversampling == 0)
				{
					{
						*AWGtablexy++ = to_dac(DACx_tag);//store scanx_tag value
						*AWGtablexy++ = to_dac(DACy_tag);//store scany_tag value
					}
					counterLoc++;
				}
				counter++;
				*Keytable++ = -1.0;
				*Keytable++ = -1.0;
			}
		}
		break;

	default:
		printf("Error:  No such ScanPattern");

	}
	
	leftover = NumSamplesPerCh - counter;
	for (int i=0; i<leftover; i++)
	{
			* Keytable++= locx;
			* Keytable++= locy;
	}
	

	//determine ArinaGate pattern table
	if (true) //before it was once  !LowRes && use_ARINA
	{
		int intlocx, intlocy;
		int16_t* LOCdata = (int16_t*)pvBufferKey; //pointer to search over key (location) data
		int record_position = 0;
		int state_counter = 0;//ArinaGate connected to AWG trig-out, otherwise= REPtrigger_to_output_delay_base* recent_oversampling; // after the trigger-in signal the output is still delayed by REPtrigger_to_output_delay_base AWG clocks
		bool record_state = false;
		bool previous_record_state = false;
		int packet_id = 0;
		int half_oversampling = (int)(recent_oversampling/2);
		for (int j = 0; j < NumSamplesPerCh; j++)
		{
			locx = *LOCdata++;
			locy = *LOCdata++;
			if (j % recent_oversampling == half_oversampling)
			{
				intlocx = (int)(locx);
				intlocy = (int)(locy);
				record_state = (intlocx >= 0 && intlocy >= 0 && intlocx < width && intlocy < height);
				if (record_state == previous_record_state && state_counter < 32765 )
				{
					state_counter++;
				}
				else
				{
					if (record_position <= Gate_pattern_count - 4)
					{
						Gate_pattern[record_position] = (byte)(state_counter & 0b11111111);
						record_position++;
						if (previous_record_state)
						{
							Gate_pattern[record_position] = (byte)((state_counter >> 8 & 0b01111111) | 0b10000000);
						}
						else
						{
							Gate_pattern[record_position] = (byte)((state_counter >> 8 & 0b01111111));
						}
						record_position++;
					}
					state_counter = 1;
				}
				previous_record_state = record_state;
			}
		}
		pattern_number_of_active_entries=record_position/2;
		//for (int pos = record_position; pos < Gate_pattern_count; pos++)
		//{
		//	Gate_pattern[pos] = 0;
		//}
	}

}

/*
**************************************************************************
Calculate DAC of scanx/y with correction of LR circuit 
**************************************************************************
*/
/*
double correctedLOCxy(double LOCxy, double delay, int thesamplerate, int counter, double * prev_correctedLOCxy)
{
	double result;
	double previous= *prev_correctedLOCxy;
	double dt=1.0/thesamplerate;
	if (counter==0) previous=LOCxy;
	result=(LOCxy+(delay/dt)*previous)/(1+(delay/dt));
	if (result>32766) result=32766;
	if (result<-32767) result=-32767;
	*prev_correctedLOCxy=result;
	return (double)result;
	
}*/

bool acquire_image(const int width, const int height, short* pData)  //pData is the output that contains clean width*height values
{
		short* CHdata = (short*)pvBufferREC; //pointer to search over arranged input data
		int16_t* LOCdata = (int16_t*)pvBufferKey; //pointer to search over key (location) data
		short* fillcounter = new short[height * width];
		//int32* fillresults = new int32[height * width];
		float* fillresults = new float[height * width];
		int16_t locx, locy; 
		int intlocx, intlocy, intloc;
		for (int j = 0; j < height * width; j++)
		{
			fillresults[j] = 0.0;
			fillcounter[j] = 0;
		}
		//determine filling points, both for pixel calculation and for building gate pattern
		int record_position = 0;
		//int state_counter = REPtrigger_to_output_delay_base * recent_oversampling; // after the active trigger signal the output is still delayed by REPtrigger_to_output_delay_base AWG clocks
		bool record_state = false;
		for (int j = 0; j < NumSamplesPerCh; j++, CHdata++)
		{
			locx = *LOCdata++;
			locy = *LOCdata++;
			intlocx = (int)(locx);
			intlocy = (int)(locy);
			if (intlocx >= 0 && intlocy >= 0 && intlocx < width && intlocy < height)
			{
				intloc = intlocx + intlocy * width;
				fillresults[intloc] += (float)(*CHdata);
				fillcounter[intloc] += 1;
			}

		}
		int counter = 0;
		float sumval = 0;
		float result_in_pixel;//, truncated_result_in_pixel;
		for (int j = 0; j < height * width; j++)
		{
			if (fillcounter[j] > 0)
			{
				result_in_pixel = fillresults[j] / (float)fillcounter[j]; //use average (all points marked with same x,y are oversampling points
				//truncated_result_in_pixel = result_in_pixel - (double)saved_bias[ch]; //remove the measured bias in zero illumination
				//if (truncated_result_in_pixel < 0) truncated_result_in_pixel = 0;
				fillresults[j] = result_in_pixel;//truncated_result_in_pixel;
				counter++;
				sumval += result_in_pixel;
			}
		}
		float meanvalue = (sumval / (float)counter);
		float ptN = 0, ptS = 0, ptE = 0, ptW = 0;
		for (int j = 0; j < height * width; j++)
		{
			if (fillcounter[j] == 0)
			{
				//truncated_result_in_pixel = (double)meanvalue[ch] - (double)saved_bias[ch];
				//if (truncated_result_in_pixel < 0) truncated_result_in_pixel = 0;
				fillresults[j] = meanvalue;//truncated_result_in_pixel;
				if (j > width + 1 && j < (height - 1) * width - 1)
				{
					ptN = fillresults[j - width];
					ptS = fillresults[j + width];
					ptE = fillresults[j + 1];
					ptW = fillresults[j - 1];
					if (ptN > 0 && ptS > 0 && ptE > 0 && ptW > 0)
						fillresults[j] = (ptN + ptS + ptE + ptW) / 4.0;
				}

			}
		}
		
		for (int j = 0; j < height * width; j++)
		{
			//image to SerialEM is short: choose the right 16bit settings in serialEM, or choose camera- devide 16 bits by 2 
			pData[j] = (short)fillresults[(j % width) + (height - (j / width) - 1) * width]; //copy image and invert y axis since serialEM flips it 
		}
		delete[] fillcounter;
		delete[] fillresults;
		return true;
}



int activate_scan(double integration_time_us, int number_of_triggers)
{

    fd_mem = ::open("/dev/mem", O_RDWR | O_SYNC);
    if (fd_mem < 0) { perror("open /dev/mem");fprintf(stderr, "Internal communication error\n"); return 1; }

	if (!LowRes && is_tilt_series) 
	{
		fprintf(stderr, "Tiltview index: %d\n",TiltIndex);
	}
	if (!LowRes && use_ARINA) {
		run_Arina_acquire(integration_time_us,number_of_triggers);
		if (!is_tilt_series) use_ARINA=false;
		usleep(4000000);//Dectris server will not accept immediatley triggers 
	}
	bool ok = run_scan_with_gate_and_adc();
	if (!LowRes && is_tilt_series) {
		TiltIndex++;
		if (TiltIndex>=max_tilts)
		{
			is_tilt_series=false;
			use_ARINA=false;
			TiltIndex=0;
		}
	}
	
	   // Cleanup (you can keep mappings for repeated programming)
    munmap((void *)bram_mmap,  SPAN_MEM);
    munmap((void *)map_ctrl_mmap, SPAN_CTRL);
	munmap((void *)ddr_ptr_mmap, dmm_size);  //Should release the volatile memory used for both ADC and DAC
	munmap((void *)dma_ctrl_mmap, DMA_SPAN);
	close(fd_mem);	

	if (!ok) {
		fprintf(stderr, "Scan failed\n");
		return 1;}

	return 0;	

}
