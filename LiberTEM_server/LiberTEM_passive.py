# LiberTEM_passive.py: A LiberTEM-live server for SavvyScan-ver-2b, written by Shahar Seifer (2023), Elbaum lab, Weizmann Institute of Science
# citation for LiberTEM-live code:  10.5281/zenodo.4916315
# citation for SSB implementation: M&M (2021), 27, 1078–1092, doi:10.1017/S1431927621012423
# citation for ptychography code: 10.5281/zenodo.5055126
# citation for SavvyScan system: M&M (2021) https://doi.org/10.1017/S1431927621012861
# citation for SavvyScan code: https://github.com/Pr4Et/SavvyScan
# part of the code obtained from: https://github.com/LiberTEM/LiberTEM-live/blob/eb12fdbfd66984c8389a0e8c4757074748af2de0/docs/source/detectors/dectris.rst
# part obtained from: https://ptychography-4-0.github.io/ptychography/algorithms/live-ssb.html#SSB-setup
# 
from libertem.viz.bqp import BQLive2DPlot
from libertem.viz.mpl import MPLLive2DPlot
from libertem_live.api import LiveContext
from libertem.udf import UDF
from libertem.corrections.coordinates import flip_y, rotate_deg, identity
from libertem.analysis import com as com_analysis
from libertem.udf.masks import ApplyMasksUDF
from libertem.common.container import MaskContainer
from ptychography40.reconstruction.ssb import SSB_UDF, generate_masks
from ptychography40.reconstruction.common import wavelength, get_shifted
from scanf import scanf
import base64
import zmq
from PIL import Image
import numpy as np

class DemoNavUDF(UDF):
    def get_result_buffers(self):
        return {'nav_sum': self.buffer(kind='nav', dtype='float32'),'nav_max':self.buffer(kind='nav', dtype='float32'),}
    def process_frame(self, frame):
        self.results.nav_sum[:] = np.sum(frame)
        self.results.nav_max[:] = np.max(frame)

def center_shifts(udf_result):
    '''
    Derive center of mass results from the UDF results
    and apply coordinate correction.
    '''
    y_centers_raw, x_centers_raw = com_analysis.center_shifts(img_sum=udf_result['intensity'].data[..., 0],img_y=udf_result['intensity'].data[..., 1],img_x=udf_result['intensity'].data[..., 2],ref_y=48,ref_x=48,)
    shape = y_centers_raw.shape
    y_centers, x_centers = rotate_deg(0) @ flip_y() @ (y_centers_raw.reshape(-1), x_centers_raw.reshape(-1))
    y_centers = y_centers.reshape(shape)
    x_centers = x_centers.reshape(shape)
    return (y_centers, x_centers)

def magnitude(udf_result, damage):
    return (com_analysis.magnitude(*center_shifts(udf_result)), damage)

#executor = #PipelinedExecutor(spec=PipelinedExecutor.make_spec(cpus=range(10), #cudas=[]))
#ctx = LiveContext(executor=executor, plot_class=BQLive2DPlot)
#ctx = LiveContext(plot_class=MPLLive2DPlot)
ctx = LiveContext(plot_class=BQLive2DPlot)
udfs=[DemoNavUDF()]

# prepare for acquisition, setting up scan parameters etc.
conn = ctx.make_connection('dectris').open(
    api_host="192.168.100.70",
    api_port=80,
    data_host="192.168.100.70",
    data_port=9999,
    buffer_size=2048,
)

context=zmq.Context()
socket=context.socket(zmq.REP) #replier
socket.bind("tcp://*:5559") #It is SavvyScan IP in local network. 

reply_txt="OK"
nx=0
ny=0
ri=0
ro=10
ptype="OFF"

with conn:
    while True:
            message=socket.recv()
            print("Received request: %s" % message)
            socket.send_string(reply_txt)
            if message[0:4]==b'EXIT':
                 ptype="EXIT"
                 break
            elif message[0:3]==b'SUM':
                 ptype="SUM"
            elif message[0:2]==b'BF':
                 ptype="BF"
                 [ri,ro]=scanf("BF;ri=%d,ro=%d",message.decode('utf-8'))
            elif message[0:3]==b'SSB':
                 #Use appropriate scan size limited by GPU memory and 1000 frames/sec processing rate
                 ptype="SSB"
                 [KV,pix_nm,scon_mrad,disc_px]=scanf("SSB;KV=%d,pix_nm=%f,scon_mrad=%f,disc_px=%d",message.decode('utf-8'))
                 rec_params = {
                    # NumPy base dtype for the calculation. Set to numpy.float64 for higher precision
                    'dtype': np.float32,
                    # Wavelength of the illuminating radiation in m (fill in KV)
                    'lamb': wavelength(KV),
                    # Scan step in x and y direction in m
                    'dpix': pix_nm*1e-9/100,
                    # Semiconvergence angle of the illumination in radians
                    'semiconv': scon_mrad*0.001,
                    # Diameter of the zero order disk on the detector in pixels
                    'semiconv_pix': disc_px,
                    # Affine transformation matrix for adjusting rotation and handedness between scan and detector coordinates.
                    # The scale is handled by the other parameters.
                    'transformation': rotate_deg(0) @ flip_y(),
                    # Position of the direct beam on the detecot in px
                    'cy': 48,
                    'cx': 48,
                    # Minimum size of a trotter. Increase this to suppress noise from very small trotters.
                    'cutoff': 10
                 }
                 cutoff_freq = np.float32('inf')
            elif message[0:3]==b'COM':
                 ptype="COM"
                 [disc_px]=scanf("COM;disc_pix=%d",message.decode('utf-8'))
                 masks = com_analysis.com_masks_factory(
                     detector_y=96,
                     detector_x=96,
                     cx=48,
                     cy=48,
                     r=disc_px/2+4,
                 )
                 com_udf = ApplyMasksUDF(masks)
            elif message[0:5]==b'START':
                [nx,ny]=scanf("START %d,%d",message.decode('utf-8'))
                pending_aq = conn.wait_for_acquisition(timeout=60)
                print("Handling acquisition")
                img = Image.fromarray(np.uint16([2,2])) 
                img.save("/home/ftpuser/LiberTEM.png") #first replace result with empty file
                aq=ctx.make_acquisition(conn=conn,nav_shape=(ny,nx),frames_per_partition=nx*ny/32,pending_aq=pending_aq,)
                if ptype=="SUM":
                    BF_RES=ctx.run_udf(dataset=aq, udf=udfs,plots=[['nav_sum'],])
                    BFw = BF_RES[0]
                    BF_VAL = BFw['nav_sum'].data
                    message=b'OFF\0'
                elif ptype=="BF":
                    ring = ctx.create_ring_analysis(dataset=aq, ri=ri, ro=ro)
                    BF_RES = ctx.run(ring, progress=True)
                    BF_VAL = BF_RES['intensity']
                    message=b'OFF\0'
                elif ptype=="SSB": 
                    mask_params = {
                        # Shape of the reconstructed area
                        'reconstruct_shape': tuple(aq.shape.nav),
                        # Shape of a detector frame
                        'mask_shape': tuple(aq.shape.sig),
                        # Use the faster shifting method to generate trotters
                        'method': 'shift',
                    }
                    trotters = generate_masks(**rec_params, **mask_params)
                    mask_container = MaskContainer(mask_factories=lambda: trotters, dtype=trotters.dtype, count=trotters.shape[0])
                    ssb_udf = SSB_UDF(**rec_params, mask_container=mask_container)
                    BF = ctx.run_udf(dataset=aq, udf=ssb_udf)
                    #plot = MPLLive2DPlot(aq, ssb_udf, channel=phase) #amplitude/ phase 
                    BF_VAL = BF['phase'].data
                    message=b'OFF\0'
                elif ptype=="COM":
                    BF = ctx.run_udf(dataset=aq, udf=com_udf)
                    BF_RES = BF['intensity']
                    #BF_RES = BQLive2DPlot(aq, com_udf, channel=magnitude) 
                    #BF_RES = MPLLive2DPlot(aq, com_udf, channel=magnitude)
                    BF_VAL1 = BF_RES.data[:,:,1]
                    BF_VAL2 = BF_RES.data[:,:,2]
                    BF_VAL = np.sqrt(np.power(BF_VAL1,2)+np.power(BF_VAL2,2))
                    message=b'OFF\0'
                elif ptype=="SSB" and nx>1024:
                    printf("Not enough resources")
                    BF_VAL = np.zeros(np.uint16([4,4]))
                    message=b'OFF\0'
                else:
                    continue
                #BF_VAL.shape=(ny,nx)
                BF_norm=(BF_VAL-np.min(BF_VAL))/(0.001+np.max(BF_VAL)-np.min(BF_VAL))
                img = Image.fromarray(np.uint16(32767 * BF_norm)) 
                img.save("/home/ftpuser/LiberTEM.png")
                print("Image saved")

conn.close()
