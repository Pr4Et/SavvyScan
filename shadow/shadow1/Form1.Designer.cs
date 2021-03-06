﻿namespace shadow1
{
    public partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.SavetoFolder = new System.Windows.Forms.Button();
            this.FolderName = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.OutputScanAmp = new System.Windows.Forms.NumericUpDown();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.InputScanAmp = new System.Windows.Forms.NumericUpDown();
            this.ScanArgument = new System.Windows.Forms.NumericUpDown();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.AspectRatio = new System.Windows.Forms.NumericUpDown();
            this.END = new System.Windows.Forms.Button();
            this.timer1 = new System.Windows.Forms.Timer(this.components);
            this.ScanMode_comboBox1 = new System.Windows.Forms.ComboBox();
            this.ServerAddress = new System.Windows.Forms.TextBox();
            this.label6 = new System.Windows.Forms.Label();
            this.pictureBox1 = new System.Windows.Forms.PictureBox();
            this.label7 = new System.Windows.Forms.Label();
            this.StepAngle = new System.Windows.Forms.NumericUpDown();
            this.TiltAngle = new System.Windows.Forms.NumericUpDown();
            this.tomogramIndex = new System.Windows.Forms.NumericUpDown();
            this.label8 = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.label10 = new System.Windows.Forms.Label();
            this.startTomo = new System.Windows.Forms.Button();
            this.NextFrame = new System.Windows.Forms.Button();
            this.EndTomo = new System.Windows.Forms.Button();
            this.chosenCH = new System.Windows.Forms.NumericUpDown();
            this.label11 = new System.Windows.Forms.Label();
            this.MonitorFile = new System.Windows.Forms.Button();
            this.Folder2Monitor = new System.Windows.Forms.TextBox();
            this.label12 = new System.Windows.Forms.Label();
            this.magnification_combo = new System.Windows.Forms.ComboBox();
            this.Defocus_view = new System.Windows.Forms.TextBox();
            this.label13 = new System.Windows.Forms.Label();
            this.checkBox1 = new System.Windows.Forms.CheckBox();
            this.saveMAT = new System.Windows.Forms.CheckBox();
            this.Save_MATKEY = new System.Windows.Forms.CheckBox();
            this.GlanceFiles = new System.Windows.Forms.Button();
            this.glancefiles_folder = new System.Windows.Forms.TextBox();
            this.glance_delay = new System.Windows.Forms.NumericUpDown();
            this.label14 = new System.Windows.Forms.Label();
            this.StartMulti = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.OutputScanAmp)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.InputScanAmp)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.ScanArgument)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.AspectRatio)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.StepAngle)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.TiltAngle)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.tomogramIndex)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.chosenCH)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.glance_delay)).BeginInit();
            this.SuspendLayout();
            // 
            // textBox1
            // 
            this.textBox1.Location = new System.Drawing.Point(641, 71);
            this.textBox1.Margin = new System.Windows.Forms.Padding(2);
            this.textBox1.Multiline = true;
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(352, 331);
            this.textBox1.TabIndex = 0;
            // 
            // SavetoFolder
            // 
            this.SavetoFolder.Location = new System.Drawing.Point(131, 32);
            this.SavetoFolder.Margin = new System.Windows.Forms.Padding(2);
            this.SavetoFolder.Name = "SavetoFolder";
            this.SavetoFolder.Size = new System.Drawing.Size(81, 36);
            this.SavetoFolder.TabIndex = 1;
            this.SavetoFolder.Text = "Save to folder:";
            this.SavetoFolder.UseVisualStyleBackColor = true;
            this.SavetoFolder.Click += new System.EventHandler(this.SavetoFolder_Click);
            // 
            // FolderName
            // 
            this.FolderName.Location = new System.Drawing.Point(15, 72);
            this.FolderName.Margin = new System.Windows.Forms.Padding(2);
            this.FolderName.Name = "FolderName";
            this.FolderName.Size = new System.Drawing.Size(343, 20);
            this.FolderName.TabIndex = 2;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(19, 122);
            this.label1.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(62, 13);
            this.label1.TabIndex = 4;
            this.label1.Text = "ScanMode:";
            // 
            // OutputScanAmp
            // 
            this.OutputScanAmp.Increment = new decimal(new int[] {
            100,
            0,
            0,
            0});
            this.OutputScanAmp.Location = new System.Drawing.Point(29, 305);
            this.OutputScanAmp.Margin = new System.Windows.Forms.Padding(2);
            this.OutputScanAmp.Maximum = new decimal(new int[] {
            5000,
            0,
            0,
            0});
            this.OutputScanAmp.Name = "OutputScanAmp";
            this.OutputScanAmp.Size = new System.Drawing.Size(118, 20);
            this.OutputScanAmp.TabIndex = 5;
            this.OutputScanAmp.Value = new decimal(new int[] {
            2000,
            0,
            0,
            0});
            this.OutputScanAmp.ValueChanged += new System.EventHandler(this.OutputScanAmp_ValueChanged);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(26, 282);
            this.label2.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(143, 13);
            this.label2.TabIndex = 6;
            this.label2.Text = "Output Scan Amplitude [mV]:";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(200, 282);
            this.label3.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(142, 13);
            this.label3.TabIndex = 7;
            this.label3.Text = "Input Range Amplitude [mV]:";
            this.label3.Click += new System.EventHandler(this.label3_Click);
            // 
            // InputScanAmp
            // 
            this.InputScanAmp.Increment = new decimal(new int[] {
            100,
            0,
            0,
            0});
            this.InputScanAmp.Location = new System.Drawing.Point(203, 305);
            this.InputScanAmp.Margin = new System.Windows.Forms.Padding(2);
            this.InputScanAmp.Maximum = new decimal(new int[] {
            5000,
            0,
            0,
            0});
            this.InputScanAmp.Name = "InputScanAmp";
            this.InputScanAmp.Size = new System.Drawing.Size(118, 20);
            this.InputScanAmp.TabIndex = 8;
            this.InputScanAmp.Value = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            this.InputScanAmp.ValueChanged += new System.EventHandler(this.InputScanAmp_ValueChanged);
            // 
            // ScanArgument
            // 
            this.ScanArgument.DecimalPlaces = 4;
            this.ScanArgument.Increment = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            this.ScanArgument.Location = new System.Drawing.Point(22, 225);
            this.ScanArgument.Margin = new System.Windows.Forms.Padding(2);
            this.ScanArgument.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.ScanArgument.Name = "ScanArgument";
            this.ScanArgument.Size = new System.Drawing.Size(118, 20);
            this.ScanArgument.TabIndex = 12;
            this.ScanArgument.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.ScanArgument.ValueChanged += new System.EventHandler(this.ScanArgument_ValueChanged);
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(19, 202);
            this.label4.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(83, 13);
            this.label4.TabIndex = 11;
            this.label4.Text = "Scan Argument:";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(200, 202);
            this.label5.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(68, 13);
            this.label5.TabIndex = 10;
            this.label5.Text = "AspectRatio:";
            // 
            // AspectRatio
            // 
            this.AspectRatio.DecimalPlaces = 4;
            this.AspectRatio.Increment = new decimal(new int[] {
            1,
            0,
            0,
            65536});
            this.AspectRatio.Location = new System.Drawing.Point(203, 225);
            this.AspectRatio.Margin = new System.Windows.Forms.Padding(2);
            this.AspectRatio.Maximum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.AspectRatio.Name = "AspectRatio";
            this.AspectRatio.Size = new System.Drawing.Size(118, 20);
            this.AspectRatio.TabIndex = 9;
            this.AspectRatio.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.AspectRatio.ValueChanged += new System.EventHandler(this.AspectRatio_ValueChanged);
            // 
            // END
            // 
            this.END.Location = new System.Drawing.Point(30, 360);
            this.END.Margin = new System.Windows.Forms.Padding(2);
            this.END.Name = "END";
            this.END.Size = new System.Drawing.Size(65, 33);
            this.END.TabIndex = 13;
            this.END.Text = "END";
            this.END.UseVisualStyleBackColor = true;
            this.END.Click += new System.EventHandler(this.END_Click);
            // 
            // timer1
            // 
            this.timer1.Enabled = true;
            this.timer1.Interval = 1000;
            this.timer1.Tick += new System.EventHandler(this.timer1_Tick);
            // 
            // ScanMode_comboBox1
            // 
            this.ScanMode_comboBox1.DisplayMember = "0";
            this.ScanMode_comboBox1.FormattingEnabled = true;
            this.ScanMode_comboBox1.Items.AddRange(new object[] {
            "Standard (zigzag)",
            "Spiral",
            "Linear Mandela: elliptic"});
            this.ScanMode_comboBox1.Location = new System.Drawing.Point(22, 144);
            this.ScanMode_comboBox1.Margin = new System.Windows.Forms.Padding(2);
            this.ScanMode_comboBox1.Name = "ScanMode_comboBox1";
            this.ScanMode_comboBox1.Size = new System.Drawing.Size(191, 21);
            this.ScanMode_comboBox1.TabIndex = 14;
            this.ScanMode_comboBox1.Text = "Standard (zigzag)";
            this.ScanMode_comboBox1.ValueMember = "0";
            this.ScanMode_comboBox1.SelectedIndexChanged += new System.EventHandler(this.ScanMode_comboBox1_SelectedIndexChanged);
            // 
            // ServerAddress
            // 
            this.ServerAddress.Enabled = false;
            this.ServerAddress.Location = new System.Drawing.Point(113, 433);
            this.ServerAddress.Margin = new System.Windows.Forms.Padding(2);
            this.ServerAddress.Name = "ServerAddress";
            this.ServerAddress.Size = new System.Drawing.Size(155, 20);
            this.ServerAddress.TabIndex = 15;
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(27, 435);
            this.label6.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(81, 13);
            this.label6.TabIndex = 16;
            this.label6.Text = "Server address:";
            // 
            // pictureBox1
            // 
            this.pictureBox1.BackColor = System.Drawing.SystemColors.GradientInactiveCaption;
            this.pictureBox1.Location = new System.Drawing.Point(434, 42);
            this.pictureBox1.Margin = new System.Windows.Forms.Padding(2);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(167, 162);
            this.pictureBox1.TabIndex = 17;
            this.pictureBox1.TabStop = false;
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(444, 24);
            this.label7.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(76, 13);
            this.label7.TabIndex = 18;
            this.label7.Text = "Beam position:";
            // 
            // StepAngle
            // 
            this.StepAngle.Location = new System.Drawing.Point(433, 320);
            this.StepAngle.Margin = new System.Windows.Forms.Padding(2);
            this.StepAngle.Name = "StepAngle";
            this.StepAngle.Size = new System.Drawing.Size(80, 20);
            this.StepAngle.TabIndex = 19;
            this.StepAngle.Value = new decimal(new int[] {
            2,
            0,
            0,
            0});
            // 
            // TiltAngle
            // 
            this.TiltAngle.Location = new System.Drawing.Point(433, 360);
            this.TiltAngle.Margin = new System.Windows.Forms.Padding(2);
            this.TiltAngle.Maximum = new decimal(new int[] {
            90,
            0,
            0,
            0});
            this.TiltAngle.Minimum = new decimal(new int[] {
            90,
            0,
            0,
            -2147483648});
            this.TiltAngle.Name = "TiltAngle";
            this.TiltAngle.Size = new System.Drawing.Size(80, 20);
            this.TiltAngle.TabIndex = 20;
            this.TiltAngle.ValueChanged += new System.EventHandler(this.TiltAngle_ValueChanged);
            // 
            // tomogramIndex
            // 
            this.tomogramIndex.Location = new System.Drawing.Point(433, 278);
            this.tomogramIndex.Margin = new System.Windows.Forms.Padding(2);
            this.tomogramIndex.Name = "tomogramIndex";
            this.tomogramIndex.Size = new System.Drawing.Size(80, 20);
            this.tomogramIndex.TabIndex = 21;
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(431, 263);
            this.label8.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(65, 13);
            this.label8.TabIndex = 22;
            this.label8.Text = "Frame Index";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(431, 305);
            this.label9.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(89, 13);
            this.label9.TabIndex = 23;
            this.label9.Text = "Angle Increments";
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(431, 345);
            this.label10.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(51, 13);
            this.label10.TabIndex = 24;
            this.label10.Text = "Tilt Angle";
            // 
            // startTomo
            // 
            this.startTomo.Location = new System.Drawing.Point(549, 270);
            this.startTomo.Margin = new System.Windows.Forms.Padding(2);
            this.startTomo.Name = "startTomo";
            this.startTomo.Size = new System.Drawing.Size(79, 31);
            this.startTomo.TabIndex = 25;
            this.startTomo.Text = "Start Tomo";
            this.startTomo.UseVisualStyleBackColor = true;
            this.startTomo.Click += new System.EventHandler(this.startTomo_Click);
            // 
            // NextFrame
            // 
            this.NextFrame.Location = new System.Drawing.Point(549, 320);
            this.NextFrame.Margin = new System.Windows.Forms.Padding(2);
            this.NextFrame.Name = "NextFrame";
            this.NextFrame.Size = new System.Drawing.Size(79, 30);
            this.NextFrame.TabIndex = 26;
            this.NextFrame.Text = "Next Frame";
            this.NextFrame.UseVisualStyleBackColor = true;
            this.NextFrame.Click += new System.EventHandler(this.NextFrame_Click);
            // 
            // EndTomo
            // 
            this.EndTomo.Location = new System.Drawing.Point(549, 369);
            this.EndTomo.Margin = new System.Windows.Forms.Padding(2);
            this.EndTomo.Name = "EndTomo";
            this.EndTomo.Size = new System.Drawing.Size(77, 33);
            this.EndTomo.TabIndex = 27;
            this.EndTomo.Text = "End Tomo";
            this.EndTomo.UseVisualStyleBackColor = true;
            this.EndTomo.Click += new System.EventHandler(this.EndTomo_Click);
            // 
            // chosenCH
            // 
            this.chosenCH.Location = new System.Drawing.Point(290, 144);
            this.chosenCH.Maximum = new decimal(new int[] {
            7,
            0,
            0,
            0});
            this.chosenCH.Name = "chosenCH";
            this.chosenCH.Size = new System.Drawing.Size(120, 20);
            this.chosenCH.TabIndex = 28;
            this.chosenCH.ValueChanged += new System.EventHandler(this.chosenCH_ValueChanged);
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(287, 122);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(126, 13);
            this.label11.TabIndex = 29;
            this.label11.Text = "Channel sent to SerialEM";
            // 
            // MonitorFile
            // 
            this.MonitorFile.Location = new System.Drawing.Point(641, 433);
            this.MonitorFile.Name = "MonitorFile";
            this.MonitorFile.Size = new System.Drawing.Size(168, 29);
            this.MonitorFile.TabIndex = 30;
            this.MonitorFile.Text = "Monitor commands via files";
            this.MonitorFile.UseVisualStyleBackColor = true;
            this.MonitorFile.Click += new System.EventHandler(this.button1_Click);
            // 
            // Folder2Monitor
            // 
            this.Folder2Monitor.Location = new System.Drawing.Point(461, 438);
            this.Folder2Monitor.Name = "Folder2Monitor";
            this.Folder2Monitor.Size = new System.Drawing.Size(165, 20);
            this.Folder2Monitor.TabIndex = 31;
            this.Folder2Monitor.Text = "d:\\shared";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(169, 363);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(80, 13);
            this.label12.TabIndex = 32;
            this.label12.Text = "Magnification X";
            // 
            // magnification_combo
            // 
            this.magnification_combo.DisplayMember = "0";
            this.magnification_combo.FormattingEnabled = true;
            this.magnification_combo.Items.AddRange(new object[] {
            "9900",
            "14000",
            "20000",
            "28000",
            "40000",
            "56000",
            "79000",
            "110000"});
            this.magnification_combo.Location = new System.Drawing.Point(254, 360);
            this.magnification_combo.Margin = new System.Windows.Forms.Padding(2);
            this.magnification_combo.Name = "magnification_combo";
            this.magnification_combo.Size = new System.Drawing.Size(120, 21);
            this.magnification_combo.TabIndex = 33;
            this.magnification_combo.Text = "9900";
            this.magnification_combo.ValueMember = "0";
            this.magnification_combo.SelectedIndexChanged += new System.EventHandler(this.magnification_combo_SelectedIndexChanged);
            // 
            // Defocus_view
            // 
            this.Defocus_view.Location = new System.Drawing.Point(254, 397);
            this.Defocus_view.Name = "Defocus_view";
            this.Defocus_view.Size = new System.Drawing.Size(129, 20);
            this.Defocus_view.TabIndex = 34;
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(198, 400);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(50, 13);
            this.label13.TabIndex = 35;
            this.label13.Text = "Defocus:";
            // 
            // checkBox1
            // 
            this.checkBox1.AutoSize = true;
            this.checkBox1.Location = new System.Drawing.Point(352, 438);
            this.checkBox1.Name = "checkBox1";
            this.checkBox1.Size = new System.Drawing.Size(63, 17);
            this.checkBox1.TabIndex = 36;
            this.checkBox1.Text = "Lothar\'s";
            this.checkBox1.UseVisualStyleBackColor = true;
            this.checkBox1.CheckedChanged += new System.EventHandler(this.checkBox1_CheckedChanged);
            // 
            // saveMAT
            // 
            this.saveMAT.AutoSize = true;
            this.saveMAT.Location = new System.Drawing.Point(295, 42);
            this.saveMAT.Name = "saveMAT";
            this.saveMAT.Size = new System.Drawing.Size(77, 17);
            this.saveMAT.TabIndex = 37;
            this.saveMAT.Text = "Save MAT";
            this.saveMAT.UseVisualStyleBackColor = true;
            this.saveMAT.CheckedChanged += new System.EventHandler(this.saveMAT_CheckedChanged);
            // 
            // Save_MATKEY
            // 
            this.Save_MATKEY.AutoSize = true;
            this.Save_MATKEY.Location = new System.Drawing.Point(378, 42);
            this.Save_MATKEY.Name = "Save_MATKEY";
            this.Save_MATKEY.Size = new System.Drawing.Size(44, 17);
            this.Save_MATKEY.TabIndex = 38;
            this.Save_MATKEY.Text = "Key";
            this.Save_MATKEY.UseVisualStyleBackColor = true;
            this.Save_MATKEY.CheckedChanged += new System.EventHandler(this.Save_MATKEY_CheckedChanged);
            // 
            // GlanceFiles
            // 
            this.GlanceFiles.Location = new System.Drawing.Point(641, 12);
            this.GlanceFiles.Name = "GlanceFiles";
            this.GlanceFiles.Size = new System.Drawing.Size(114, 37);
            this.GlanceFiles.TabIndex = 39;
            this.GlanceFiles.Text = "Glance through files";
            this.GlanceFiles.UseVisualStyleBackColor = true;
            this.GlanceFiles.Click += new System.EventHandler(this.GlanceFiles_Click);
            // 
            // glancefiles_folder
            // 
            this.glancefiles_folder.Location = new System.Drawing.Point(760, 29);
            this.glancefiles_folder.Margin = new System.Windows.Forms.Padding(2);
            this.glancefiles_folder.Name = "glancefiles_folder";
            this.glancefiles_folder.Size = new System.Drawing.Size(233, 20);
            this.glancefiles_folder.TabIndex = 40;
            this.glancefiles_folder.TextChanged += new System.EventHandler(this.glancefiles_folder_TextChanged);
            // 
            // glance_delay
            // 
            this.glance_delay.Increment = new decimal(new int[] {
            100,
            0,
            0,
            0});
            this.glance_delay.Location = new System.Drawing.Point(920, 4);
            this.glance_delay.Maximum = new decimal(new int[] {
            5000,
            0,
            0,
            0});
            this.glance_delay.Name = "glance_delay";
            this.glance_delay.Size = new System.Drawing.Size(72, 20);
            this.glance_delay.TabIndex = 41;
            this.glance_delay.Value = new decimal(new int[] {
            500,
            0,
            0,
            0});
            this.glance_delay.ValueChanged += new System.EventHandler(this.glance_delay_ValueChanged);
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.Location = new System.Drawing.Point(851, 6);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(63, 13);
            this.label14.TabIndex = 42;
            this.label14.Text = "Glance time";
            // 
            // StartMulti
            // 
            this.StartMulti.Location = new System.Drawing.Point(549, 225);
            this.StartMulti.Margin = new System.Windows.Forms.Padding(2);
            this.StartMulti.Name = "StartMulti";
            this.StartMulti.Size = new System.Drawing.Size(79, 31);
            this.StartMulti.TabIndex = 43;
            this.StartMulti.Text = "Start Multi";
            this.StartMulti.UseVisualStyleBackColor = true;
            this.StartMulti.Click += new System.EventHandler(this.StartMulti_Click_1);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1004, 514);
            this.Controls.Add(this.StartMulti);
            this.Controls.Add(this.label14);
            this.Controls.Add(this.glance_delay);
            this.Controls.Add(this.glancefiles_folder);
            this.Controls.Add(this.GlanceFiles);
            this.Controls.Add(this.Save_MATKEY);
            this.Controls.Add(this.saveMAT);
            this.Controls.Add(this.checkBox1);
            this.Controls.Add(this.label13);
            this.Controls.Add(this.Defocus_view);
            this.Controls.Add(this.magnification_combo);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.Folder2Monitor);
            this.Controls.Add(this.MonitorFile);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.chosenCH);
            this.Controls.Add(this.EndTomo);
            this.Controls.Add(this.NextFrame);
            this.Controls.Add(this.startTomo);
            this.Controls.Add(this.label10);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.tomogramIndex);
            this.Controls.Add(this.TiltAngle);
            this.Controls.Add(this.StepAngle);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.pictureBox1);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.ServerAddress);
            this.Controls.Add(this.ScanMode_comboBox1);
            this.Controls.Add(this.END);
            this.Controls.Add(this.ScanArgument);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.AspectRatio);
            this.Controls.Add(this.InputScanAmp);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.OutputScanAmp);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.FolderName);
            this.Controls.Add(this.SavetoFolder);
            this.Controls.Add(this.textBox1);
            this.Margin = new System.Windows.Forms.Padding(2);
            this.Name = "Form1";
            this.Text = "Form1";
            ((System.ComponentModel.ISupportInitialize)(this.OutputScanAmp)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.InputScanAmp)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.ScanArgument)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.AspectRatio)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.StepAngle)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.TiltAngle)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.tomogramIndex)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.chosenCH)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.glance_delay)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        public System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.Button SavetoFolder;
        private System.Windows.Forms.TextBox FolderName;
        private System.Windows.Forms.Label label1;
        public System.Windows.Forms.NumericUpDown OutputScanAmp;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        public System.Windows.Forms.NumericUpDown InputScanAmp;
        private System.Windows.Forms.NumericUpDown ScanArgument;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.NumericUpDown AspectRatio;
        private System.Windows.Forms.Button END;
        private System.Windows.Forms.Timer timer1;
        public System.Windows.Forms.ComboBox ScanMode_comboBox1;
        public System.Windows.Forms.TextBox ServerAddress;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.PictureBox pictureBox1;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.NumericUpDown StepAngle;
        private System.Windows.Forms.NumericUpDown TiltAngle;
        private System.Windows.Forms.NumericUpDown tomogramIndex;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.Button startTomo;
        private System.Windows.Forms.Button NextFrame;
        private System.Windows.Forms.Button EndTomo;
        private System.Windows.Forms.NumericUpDown chosenCH;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Button MonitorFile;
        private System.Windows.Forms.TextBox Folder2Monitor;
        private System.Windows.Forms.Label label12;
        public System.Windows.Forms.ComboBox magnification_combo;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.CheckBox checkBox1;
        public System.Windows.Forms.TextBox Defocus_view;
        private System.Windows.Forms.CheckBox saveMAT;
        private System.Windows.Forms.CheckBox Save_MATKEY;
        private System.Windows.Forms.Button GlanceFiles;
        private System.Windows.Forms.TextBox glancefiles_folder;
        private System.Windows.Forms.NumericUpDown glance_delay;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.Button StartMulti;
    }
}

