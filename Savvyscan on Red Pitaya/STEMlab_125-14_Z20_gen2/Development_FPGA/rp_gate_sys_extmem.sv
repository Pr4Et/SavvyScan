`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
//
// Engineer: Shahar Seifer & Microsoft Copilot
//
//////////////////////////////////////////////////////////////////////////////////


module rp_gate_sys_extmem #(
  parameter integer ENTRY_AW = 15
)(
  // -------------------------
  // sysbus
  // -------------------------
  input  wire              clk,
  input  wire              rstn,
  input  wire [19:0]       sys_addr,
  input  wire [31:0]       sys_wdata,
  input  wire              sys_wen,
  input  wire              sys_ren,
  output reg  [31:0]       sys_rdata,
  output reg               sys_err,
  output reg               sys_ack,

  // -------------------------
  // main scan clock (must match DMA m_axis_aclk)
  // -------------------------
  input  wire              adc_clk,

  // -------------------------
  // Pattern table interface
  // -------------------------
  output reg [ENTRY_AW-1:0] pattern_idx_r,
  input  wire [15:0]        pattern_rdata,

  // -------------------------
  // DIO pulse outputs
  // -------------------------
  output reg               rp_dio0_p,
  output reg               rp_dio3_p,

  // ======================================================
  // AXI-Stream input from DMA (MM2S)
  // ======================================================
  input  wire [31:0]       s_axis_xy_tdata,
  input  wire              s_axis_xy_tvalid,
  input  wire              s_axis_xy_tlast,
  output wire              s_axis_xy_tready,

  // -------------------------
  // DAC output
  // -------------------------
  output wire signed [15:0] dac_x,
  output wire signed [15:0] dac_y
);

  // ============================
  // Register map
  // ============================
  localparam CTRL_ADDR       = 20'h00000;
  localparam TICKDIV_ADDR    = 20'h00004;
  localparam PWIDTH_ADDR     = 20'h00008;
  localparam PLEN_ADDR       = 20'h0000C;
  localparam MASK_COUNT_ADDR = 20'h00010;
  localparam MASK_STATE_ADDR = 20'h00014;

  localparam PARK_EN_ADDR    = 20'h0001C;
  localparam PARK_X_ADDR     = 20'h00020;
  localparam PARK_Y_ADDR     = 20'h00024;

  localparam XLEN_ADDR       = 20'h00028;   // # of XY pairs
  localparam TRIGGER_DELAY_ADDR     = 20'h0002C;
  
  reg        reg_arm, reg_stop, reg_send_camera;
  reg [31:0] reg_tick_div, reg_pwidth;
  reg [31:0] reg_trigger_delay;
  reg [15:0] reg_plen, reg_mask_count, reg_mask_state;
  reg        park_en;
  reg signed [15:0] park_x, park_y;
  reg [31:0] reg_xlen_pairs;
  
  localparam S_STARTUP = 1'd0;
  localparam S_NORMAL  = 1'd1;
  //localparam S_DONE    = 2'd2;
  reg axis_state;

  // =========================================================
  // sysbus read/write
  // =========================================================
  always @(posedge clk) begin
    if (!rstn) begin
      reg_arm         <= 1'b0;
      reg_stop        <= 1'b0;
      reg_send_camera <= 1'b1;
      reg_tick_div    <= 32'd1;
      reg_pwidth      <= 32'd375;
      reg_trigger_delay <= 32'd1;
      reg_plen        <= 16'd0;
      reg_mask_count  <= 16'h7FFF;
      reg_mask_state  <= 16'h8000;

      park_en         <= 1'b1;
      park_x          <= 16'sd0;
      park_y          <= 16'sd0;

      reg_xlen_pairs  <= 32'd0;

      sys_ack         <= 1'b0;
      sys_err         <= 1'b0;
      sys_rdata       <= 32'd0;
    end else begin
      sys_ack <= 1'b0;
      sys_err <= 1'b0;

      if (sys_wen) begin
        case (sys_addr)
          CTRL_ADDR:       {reg_send_camera, reg_stop, reg_arm}
                            <= sys_wdata[2:0];
          TICKDIV_ADDR:    reg_tick_div    <= sys_wdata;
          PWIDTH_ADDR:     reg_pwidth      <= sys_wdata;
          TRIGGER_DELAY_ADDR:reg_trigger_delay <= (sys_wdata < 32'd1) ? 32'd1 : sys_wdata;
          PLEN_ADDR:       reg_plen        <= sys_wdata[15:0];
          MASK_COUNT_ADDR: reg_mask_count  <= sys_wdata[15:0];
          MASK_STATE_ADDR: reg_mask_state  <= sys_wdata[15:0];

          PARK_EN_ADDR:    park_en         <= sys_wdata[0];
          PARK_X_ADDR:     park_x          <= sys_wdata[15:0];
          PARK_Y_ADDR:     park_y          <= sys_wdata[15:0];

          XLEN_ADDR:       reg_xlen_pairs  <= sys_wdata;    // # XY pairs
        endcase

        sys_ack <= 1'b1;
      end
      else if (sys_ren) begin
        case (sys_addr)
          CTRL_ADDR:       sys_rdata <= {29'd0, reg_send_camera, reg_stop, reg_arm};
          TICKDIV_ADDR:    sys_rdata <= reg_tick_div;
          PWIDTH_ADDR:     sys_rdata <= reg_pwidth;
          TRIGGER_DELAY_ADDR:sys_rdata <= reg_trigger_delay;
          PLEN_ADDR:       sys_rdata <= {16'd0, reg_plen};
          MASK_COUNT_ADDR: sys_rdata <= {16'd0, reg_mask_count};
          MASK_STATE_ADDR: sys_rdata <= {16'd0, reg_mask_state};
          PARK_EN_ADDR:    sys_rdata <= {31'd0, park_en};
          PARK_X_ADDR:     sys_rdata <= {16'd0, park_x};
          PARK_Y_ADDR:     sys_rdata <= {16'd0, park_y};
          XLEN_ADDR:       sys_rdata <= reg_xlen_pairs;
          default:         sys_rdata <= 32'd0;
        endcase

        sys_ack <= 1'b1;
      end
    end
  end

  // =========================================================
  // CDC into adc_clk domain
  // =========================================================
  reg arm_d1, arm_d2, stop_d1, stop_d2;
  always @(posedge adc_clk) begin
    arm_d1  <= reg_arm;   arm_d2  <= arm_d1;
    stop_d1 <= reg_stop;  stop_d2 <= stop_d1;
  end

  wire arm_pulse = arm_d1 & ~arm_d2;
  wire stop_sync = stop_d2;

  reg  [31:0] tick_div_q, pwidth_q, trigger_delay_q;
  reg  [15:0] plen_q, mask_cnt_q, mask_st_q;
  reg         send_cam_q;
  reg  [31:0] xlen_pairs_q;

  always @(posedge adc_clk) begin
    tick_div_q    <= reg_tick_div;
    pwidth_q      <= reg_pwidth;
    trigger_delay_q <= reg_trigger_delay;
    plen_q        <= reg_plen;
    mask_cnt_q    <= reg_mask_count;
    mask_st_q     <= reg_mask_state;
    send_cam_q    <= reg_send_camera;
    xlen_pairs_q  <= reg_xlen_pairs;
  end


  // =========================================================
  // Pattern FSM 
  // =========================================================
  reg tick_en;
  reg [31:0] tick_cnt;
  reg [14:0] countdown;
  reg        record_state, record_state_next;
  reg        pulse_active;
  reg [31:0] pulse_cnt;
  reg        first_record_after_arm;
 reg trigger_req;
 reg delay_active;
 reg [31:0] delay_cnt;
 reg pulse3_active;
 reg [31:0] pulse3_cnt;

  always @(posedge adc_clk) begin
    if (!rstn || stop_sync) begin
      tick_en       <= 1'b0;
      tick_cnt      <= 32'd0;
      countdown     <= 15'd0;
      record_state  <= 1'b0;
      axis_state <=S_STARTUP;
      rp_dio0_p     <= 1'b0;
      pulse_active  <= 1'b0;
      pulse_cnt     <= 32'd0;
      pattern_idx_r <= 0;
      first_record_after_arm <= 1'b0;
      
      trigger_req   <= 1'b0;
      delay_active  <= 1'b0;
      delay_cnt     <= 32'd0;
      pulse3_active <= 1'b0;
      pulse3_cnt    <= 32'd0;
      rp_dio3_p     <= 1'b0;

    end else begin
      if (arm_pulse) begin
        tick_en       <= 1'b1;
        tick_cnt      <= 0;
        countdown     <= 0;
        pattern_idx_r <= 0;
        record_state  <= 0;
        // Pulse 1
        rp_dio0_p     <= 1'b1;
        pulse_active  <= 1'b1;
        pulse_cnt     <= pwidth_q;

        first_record_after_arm <= 1'b1;
      end

      if (axis_state == S_STARTUP && s_axis_xy_tvalid)
            axis_state <= S_NORMAL;


      if (tick_en) begin
        if (tick_cnt == (tick_div_q - 1)) begin
          tick_cnt <= 0;

          record_state_next = record_state;
          if (countdown != 0) begin
            countdown <= countdown - 1;
          end else begin
            if (pattern_idx_r < plen_q) begin
              record_state_next = ((pattern_rdata & mask_st_q)!=0);
              pattern_idx_r    <= pattern_idx_r + 1;
              if (((pattern_rdata & mask_st_q) != 0)) begin
                    // pulse region: full length
                    countdown <= (pattern_rdata & mask_cnt_q);
              end else begin
                    // no-pulse region: subtract the tick just consumed
                    if ((pattern_rdata & mask_cnt_q) > 2)
                        countdown <= (pattern_rdata & mask_cnt_q) - 2;  //STATE_LATENCY = 2
                    else
                        countdown <= 0;
              end              
                              
            end else begin
              record_state_next = 1'b0;
            end
          end

          if (record_state_next && send_cam_q && (countdown!=0)) begin
            trigger_req <= 1'b1;
            pulse_active <= 1'b1;
            pulse_cnt    <= pwidth_q;
          end

          if (first_record_after_arm && (record_state_next && !record_state)) begin
              first_record_after_arm <= 1'b0;
          end

          record_state <= record_state_next;

        end else begin
          tick_cnt <= tick_cnt + 1;
        end
      end

      if (pulse_active) begin
        if (pulse_cnt==0) begin
          rp_dio0_p <= 1'b0;
          pulse_active <= 1'b0;
        end else begin
          pulse_cnt <= pulse_cnt - 1;
        end
      end
      
        if (trigger_req && !delay_active && !pulse3_active) begin
         delay_active <= 1'b1;
         delay_cnt <= trigger_delay_q;
         trigger_req <= 1'b0;
       end
    
       if (delay_active) begin
         if (delay_cnt <= 1) begin
           delay_active <= 1'b0;
           pulse3_active <= 1'b1;
           pulse3_cnt <= pwidth_q;
           rp_dio3_p <= 1'b1;
         end else begin
           delay_cnt <= delay_cnt - 1;
         end
       end
    
       if (pulse3_active) begin
         if (pulse3_cnt <= 1) begin
           pulse3_active <= 1'b0;
           rp_dio3_p <= 1'b0;
         end else begin
           pulse3_cnt <= pulse3_cnt - 1;
         end
       end
         
          
      
    end
  end

        
    // ==========================================
    // AXIS input staging register (1-word buffer)
    // ==========================================
    reg [31:0] pair_count;
    reg signed [15:0] x_samp, y_samp; 
    //------------------------------------------------------------
    // AXIS consumer for slow DAC with STARTUP-READY FSM
    //------------------------------------------------------------
    reg [31:0] stage0_data, stage1_data;
    reg        stage0_valid, stage1_valid;
    // Startup: READY must stay HIGH until first sample arrives
    wire startup_ready = 1'b1;
    wire normal_ready  = !stage1_valid;
    // Done: READY=1 (ignore incoming stream)
    wire done_ready    = 1'b1;
    assign s_axis_xy_tready = (axis_state == S_STARTUP) ? startup_ready :normal_ready;
    
    // "done" from DAC side
    wire done_xy = (pair_count >= xlen_pairs_q);
    reg stage0_last, stage1_last;   
    
    always @(posedge adc_clk) begin
        if (!rstn || arm_pulse) begin
            stage0_valid <= 1'b0;
            stage1_valid <= 1'b0;
            pair_count <= 0;
            stage0_last  <= 0;
            stage1_last  <= 0;
        end
        else begin
             //--------------------------------------------------------
            // AXIS capture (push)
            //--------------------------------------------------------
            if (s_axis_xy_tvalid && s_axis_xy_tready) begin
                if (!stage0_valid) begin
                    stage0_data  <= s_axis_xy_tdata;
                    stage0_valid <= 1'b1;
                    stage0_last  <= s_axis_xy_tlast;
                end
                else if (!stage1_valid) begin
                    stage1_data  <= s_axis_xy_tdata;
                    stage1_valid <= 1'b1;
                    stage1_last  <= s_axis_xy_tlast;
                end
                // else: full, backpressure prevents entry
            end
    
            //--------------------------------------------------------
            // DAC consumption (pull)
            //--------------------------------------------------------
            if (tick_en && tick_cnt == tick_div_q - 1) begin
                if (stage0_valid && !done_xy) begin
                    {y_samp, x_samp} <= stage0_data;
                    // Shift stage1 → stage0
                    stage0_data  <= stage1_data;
                    stage0_valid <= stage1_valid;
                    stage0_last      <= stage1_last;
                    stage1_valid <= 1'b0;
                    stage1_last  <= 1'b0;
                    pair_count   <= pair_count + 1;
                    if (stage0_last) begin
                        // Mark "done" immediately
                        pair_count <= xlen_pairs_q;  // Force done_xy = 1
                    end
                end
            end
        end
    end
    // =========================================================
    // DAC outputs
    // =========================================================
    assign dac_x = (park_en ) ? park_x : x_samp;
    assign dac_y = (park_en ) ? park_y : y_samp;
    
endmodule
