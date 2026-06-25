namespace shadow1
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
            this.BiasP = new System.Windows.Forms.NumericUpDown();
            this.label15 = new System.Windows.Forms.Label();
            this.CheckMag = new System.Windows.Forms.Button();
            this.magnification_select = new System.Windows.Forms.NumericUpDown();
            this.alignDif = new System.Windows.Forms.Button();
            this.zerobeam = new System.Windows.Forms.Button();
            this.label16 = new System.Windows.Forms.Label();
            this.LowResScanTime = new System.Windows.Forms.NumericUpDown();
            this.arinaON = new System.Windows.Forms.Button();
            this.PreName = new System.Windows.Forms.TextBox();
            this.Download4D = new System.Windows.Forms.Button();
            this.label19 = new System.Windows.Forms.Label();
            this.maxNumSlices = new System.Windows.Forms.NumericUpDown();
            this.label17 = new System.Windows.Forms.Label();
            this.libertem_UI = new System.Windows.Forms.ComboBox();
            this.label8 = new System.Windows.Forms.Label();
            this.send_libertem = new System.Windows.Forms.Button();
            this.chs_blocked = new System.Windows.Forms.CheckBox();
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
            ((System.ComponentModel.ISupportInitialize)(this.BiasP)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.magnification_select)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.LowResScanTime)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.maxNumSlices)).BeginInit();
            this.SuspendLayout();
            // 
            // textBox1
            // 
            this.textBox1.Location = new System.Drawing.Point(641, 71);
            this.textBox1.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.textBox1.Multiline = true;
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(352, 331);
            this.textBox1.TabIndex = 0;
            // 
            // SavetoFolder
            // 
            this.SavetoFolder.Location = new System.Drawing.Point(224, 29);
            this.SavetoFolder.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.SavetoFolder.Name = "SavetoFolder";
            this.SavetoFolder.Size = new System.Drawing.Size(81, 36);
            this.SavetoFolder.TabIndex = 1;
            this.SavetoFolder.Text = "Save to folder:";
            this.SavetoFolder.UseVisualStyleBackColor = true;
            this.SavetoFolder.Click += new System.EventHandler(this.SavetoFolder_Click);
            // 
            // FolderName
            // 
            this.FolderName.Location = new System.Drawing.Point(141, 72);
            this.FolderName.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.FolderName.Name = "FolderName";
            this.FolderName.Size = new System.Drawing.Size(216, 20);
            this.FolderName.TabIndex = 2;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(18, 112);
            this.label1.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(69, 13);
            this.label1.TabIndex = 4;
            this.label1.Text = "Scan Pattern";
            // 
            // OutputScanAmp
            // 
            this.OutputScanAmp.Increment = new decimal(new int[] {
            100,
            0,
            0,
            0});
            this.OutputScanAmp.Location = new System.Drawing.Point(29, 305);
            this.OutputScanAmp.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.InputScanAmp.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.ScanArgument.Location = new System.Drawing.Point(22, 257);
            this.ScanArgument.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.label4.Location = new System.Drawing.Point(19, 234);
            this.label4.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(83, 13);
            this.label4.TabIndex = 11;
            this.label4.Text = "Scan Argument:";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(200, 234);
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
            this.AspectRatio.Location = new System.Drawing.Point(203, 257);
            this.AspectRatio.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.END.Location = new System.Drawing.Point(30, 369);
            this.END.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            "Linear Mandela: elliptic",
            "Lissajous",
            "4X4",
            "Align Diffraction Shift"});
            this.ScanMode_comboBox1.Location = new System.Drawing.Point(21, 135);
            this.ScanMode_comboBox1.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.ServerAddress.Location = new System.Drawing.Point(113, 442);
            this.ServerAddress.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.ServerAddress.Name = "ServerAddress";
            this.ServerAddress.Size = new System.Drawing.Size(155, 20);
            this.ServerAddress.TabIndex = 15;
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(27, 444);
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
            this.pictureBox1.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(167, 162);
            this.pictureBox1.TabIndex = 17;
            this.pictureBox1.TabStop = false;
            this.pictureBox1.Click += new System.EventHandler(this.pictureBox1_Click);
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
            this.StepAngle.Location = new System.Drawing.Point(430, 369);
            this.StepAngle.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.TiltAngle.Location = new System.Drawing.Point(430, 409);
            this.TiltAngle.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.tomogramIndex.Location = new System.Drawing.Point(372, 71);
            this.tomogramIndex.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.tomogramIndex.Maximum = new decimal(new int[] {
            300,
            0,
            0,
            0});
            this.tomogramIndex.Minimum = new decimal(new int[] {
            2,
            0,
            0,
            -2147483648});
            this.tomogramIndex.Name = "tomogramIndex";
            this.tomogramIndex.Size = new System.Drawing.Size(49, 20);
            this.tomogramIndex.TabIndex = 21;
            this.tomogramIndex.ValueChanged += new System.EventHandler(this.tomogramIndex_ValueChanged);
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(428, 354);
            this.label9.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(89, 13);
            this.label9.TabIndex = 23;
            this.label9.Text = "Angle Increments";
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(428, 394);
            this.label10.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(51, 13);
            this.label10.TabIndex = 24;
            this.label10.Text = "Tilt Angle";
            // 
            // startTomo
            // 
            this.startTomo.Location = new System.Drawing.Point(521, 354);
            this.startTomo.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.startTomo.Name = "startTomo";
            this.startTomo.Size = new System.Drawing.Size(107, 31);
            this.startTomo.TabIndex = 25;
            this.startTomo.Text = "Start manual Tomo";
            this.startTomo.UseVisualStyleBackColor = true;
            this.startTomo.Click += new System.EventHandler(this.startTomo_Click);
            // 
            // NextFrame
            // 
            this.NextFrame.Location = new System.Drawing.Point(535, 396);
            this.NextFrame.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.NextFrame.Name = "NextFrame";
            this.NextFrame.Size = new System.Drawing.Size(79, 30);
            this.NextFrame.TabIndex = 26;
            this.NextFrame.Text = "Next Frame";
            this.NextFrame.UseVisualStyleBackColor = true;
            this.NextFrame.Click += new System.EventHandler(this.NextFrame_Click);
            // 
            // EndTomo
            // 
            this.EndTomo.Location = new System.Drawing.Point(547, 262);
            this.EndTomo.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.EndTomo.Name = "EndTomo";
            this.EndTomo.Size = new System.Drawing.Size(77, 33);
            this.EndTomo.TabIndex = 27;
            this.EndTomo.Text = "End Tomo";
            this.EndTomo.UseVisualStyleBackColor = true;
            this.EndTomo.Click += new System.EventHandler(this.EndTomo_Click);
            // 
            // chosenCH
            // 
            this.chosenCH.Location = new System.Drawing.Point(292, 135);
            this.chosenCH.Maximum = new decimal(new int[] {
            8,
            0,
            0,
            0});
            this.chosenCH.Name = "chosenCH";
            this.chosenCH.Size = new System.Drawing.Size(120, 20);
            this.chosenCH.TabIndex = 28;
            this.chosenCH.Value = new decimal(new int[] {
            7,
            0,
            0,
            0});
            this.chosenCH.ValueChanged += new System.EventHandler(this.chosenCH_ValueChanged);
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(289, 113);
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
            this.Folder2Monitor.Location = new System.Drawing.Point(461, 445);
            this.Folder2Monitor.Name = "Folder2Monitor";
            this.Folder2Monitor.Size = new System.Drawing.Size(165, 20);
            this.Folder2Monitor.TabIndex = 31;
            this.Folder2Monitor.Text = "d:\\shared";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(169, 372);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(80, 13);
            this.label12.TabIndex = 32;
            this.label12.Text = "Magnification X";
            // 
            // Defocus_view
            // 
            this.Defocus_view.Location = new System.Drawing.Point(254, 406);
            this.Defocus_view.Name = "Defocus_view";
            this.Defocus_view.Size = new System.Drawing.Size(129, 20);
            this.Defocus_view.TabIndex = 34;
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(198, 409);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(50, 13);
            this.label13.TabIndex = 35;
            this.label13.Text = "Defocus:";
            // 
            // checkBox1
            // 
            this.checkBox1.AutoSize = true;
            this.checkBox1.Location = new System.Drawing.Point(352, 447);
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
            this.saveMAT.Location = new System.Drawing.Point(303, 6);
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
            this.Save_MATKEY.Location = new System.Drawing.Point(380, 6);
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
            this.glancefiles_folder.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
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
            this.StartMulti.Location = new System.Drawing.Point(547, 225);
            this.StartMulti.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.StartMulti.Name = "StartMulti";
            this.StartMulti.Size = new System.Drawing.Size(79, 31);
            this.StartMulti.TabIndex = 43;
            this.StartMulti.Text = "Start Multi";
            this.StartMulti.UseVisualStyleBackColor = true;
            this.StartMulti.Click += new System.EventHandler(this.StartMulti_Click_1);
            // 
            // BiasP
            // 
            this.BiasP.Location = new System.Drawing.Point(30, 333);
            this.BiasP.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.BiasP.Maximum = new decimal(new int[] {
            50,
            0,
            0,
            0});
            this.BiasP.Minimum = new decimal(new int[] {
            50,
            0,
            0,
            -2147483648});
            this.BiasP.Name = "BiasP";
            this.BiasP.Size = new System.Drawing.Size(87, 20);
            this.BiasP.TabIndex = 44;
            this.BiasP.Value = new decimal(new int[] {
            16,
            0,
            0,
            -2147483648});
            this.BiasP.ValueChanged += new System.EventHandler(this.BiasP_ValueChanged);
            // 
            // label15
            // 
            this.label15.AutoSize = true;
            this.label15.Location = new System.Drawing.Point(121, 337);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(51, 13);
            this.label15.TabIndex = 45;
            this.label15.Text = "correct %";
            // 
            // CheckMag
            // 
            this.CheckMag.Location = new System.Drawing.Point(254, 337);
            this.CheckMag.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.CheckMag.Name = "CheckMag";
            this.CheckMag.Size = new System.Drawing.Size(120, 21);
            this.CheckMag.TabIndex = 46;
            this.CheckMag.Text = "Sync with SerialEM";
            this.CheckMag.UseVisualStyleBackColor = true;
            this.CheckMag.Click += new System.EventHandler(this.CheckMag_Click);
            // 
            // magnification_select
            // 
            this.magnification_select.Location = new System.Drawing.Point(255, 369);
            this.magnification_select.Maximum = new decimal(new int[] {
            100000000,
            0,
            0,
            0});
            this.magnification_select.Name = "magnification_select";
            this.magnification_select.Size = new System.Drawing.Size(120, 20);
            this.magnification_select.TabIndex = 47;
            this.magnification_select.Value = new decimal(new int[] {
            7000,
            0,
            0,
            0});
            this.magnification_select.ValueChanged += new System.EventHandler(this.magnification_select_ValueChanged);
            // 
            // alignDif
            // 
            this.alignDif.Location = new System.Drawing.Point(524, 18);
            this.alignDif.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.alignDif.Name = "alignDif";
            this.alignDif.Size = new System.Drawing.Size(101, 19);
            this.alignDif.TabIndex = 48;
            this.alignDif.Text = "Calibrate dif.disk";
            this.alignDif.UseVisualStyleBackColor = true;
            this.alignDif.Click += new System.EventHandler(this.alignDif_Click);
            // 
            // zerobeam
            // 
            this.zerobeam.Location = new System.Drawing.Point(14, 6);
            this.zerobeam.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.zerobeam.Name = "zerobeam";
            this.zerobeam.Size = new System.Drawing.Size(81, 43);
            this.zerobeam.TabIndex = 49;
            this.zerobeam.Text = "Zero-beam calibration";
            this.zerobeam.UseVisualStyleBackColor = true;
            this.zerobeam.Click += new System.EventHandler(this.zerobeam_Click);
            // 
            // label16
            // 
            this.label16.AutoSize = true;
            this.label16.Location = new System.Drawing.Point(461, 206);
            this.label16.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label16.Name = "label16";
            this.label16.RightToLeft = System.Windows.Forms.RightToLeft.Yes;
            this.label16.Size = new System.Drawing.Size(102, 13);
            this.label16.TabIndex = 50;
            this.label16.Text = "click to align dif.disk";
            // 
            // LowResScanTime
            // 
            this.LowResScanTime.Location = new System.Drawing.Point(371, 257);
            this.LowResScanTime.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.LowResScanTime.Name = "LowResScanTime";
            this.LowResScanTime.Size = new System.Drawing.Size(52, 20);
            this.LowResScanTime.TabIndex = 51;
            this.LowResScanTime.Value = new decimal(new int[] {
            8,
            0,
            0,
            0});
            this.LowResScanTime.ValueChanged += new System.EventHandler(this.numericUpDown1_ValueChanged);
            // 
            // arinaON
            // 
            this.arinaON.Location = new System.Drawing.Point(124, 32);
            this.arinaON.Name = "arinaON";
            this.arinaON.Size = new System.Drawing.Size(89, 29);
            this.arinaON.TabIndex = 53;
            this.arinaON.Text = "Arm 4D scan";
            this.arinaON.UseVisualStyleBackColor = true;
            this.arinaON.Click += new System.EventHandler(this.arinaON_Click_1);
            // 
            // PreName
            // 
            this.PreName.Location = new System.Drawing.Point(14, 72);
            this.PreName.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.PreName.Name = "PreName";
            this.PreName.Size = new System.Drawing.Size(125, 20);
            this.PreName.TabIndex = 54;
            // 
            // Download4D
            // 
            this.Download4D.Location = new System.Drawing.Point(337, 36);
            this.Download4D.Name = "Download4D";
            this.Download4D.Size = new System.Drawing.Size(84, 23);
            this.Download4D.TabIndex = 55;
            this.Download4D.Text = "Download 4D";
            this.Download4D.UseVisualStyleBackColor = true;
            this.Download4D.Click += new System.EventHandler(this.Download4D_Click);
            // 
            // label19
            // 
            this.label19.AutoSize = true;
            this.label19.Location = new System.Drawing.Point(356, 234);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(76, 13);
            this.label19.TabIndex = 59;
            this.label19.Text = "Thresh.time[s]:";
            // 
            // maxNumSlices
            // 
            this.maxNumSlices.Location = new System.Drawing.Point(473, 262);
            this.maxNumSlices.Maximum = new decimal(new int[] {
            121,
            0,
            0,
            0});
            this.maxNumSlices.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.maxNumSlices.Name = "maxNumSlices";
            this.maxNumSlices.Size = new System.Drawing.Size(59, 20);
            this.maxNumSlices.TabIndex = 60;
            this.maxNumSlices.Value = new decimal(new int[] {
            61,
            0,
            0,
            0});
            // 
            // label17
            // 
            this.label17.AutoSize = true;
            this.label17.Location = new System.Drawing.Point(459, 243);
            this.label17.Name = "label17";
            this.label17.Size = new System.Drawing.Size(73, 13);
            this.label17.TabIndex = 61;
            this.label17.Text = "max no. slices";
            // 
            // libertem_UI
            // 
            this.libertem_UI.DisplayMember = "0";
            this.libertem_UI.FormattingEnabled = true;
            this.libertem_UI.Items.AddRange(new object[] {
            "OFF",
            "SUM",
            "BF;ri=0,ro=10",
            "COM;disc_pix=40",
            "SSB;KV=200,pix_nm=2.1,scon_mrad=0.8,disc_px=40",
            "EXIT",
            "CONNECT"});
            this.libertem_UI.Location = new System.Drawing.Point(24, 197);
            this.libertem_UI.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.libertem_UI.Name = "libertem_UI";
            this.libertem_UI.Size = new System.Drawing.Size(191, 21);
            this.libertem_UI.TabIndex = 63;
            this.libertem_UI.Text = "OFF";
            this.libertem_UI.ValueMember = "0";
            this.libertem_UI.SelectedIndexChanged += new System.EventHandler(this.libertem_UI_SelectedIndexChanged);
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(21, 175);
            this.label8.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(56, 13);
            this.label8.TabIndex = 62;
            this.label8.Text = "LiberTEM:";
            // 
            // send_libertem
            // 
            this.send_libertem.Location = new System.Drawing.Point(224, 197);
            this.send_libertem.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.send_libertem.Name = "send_libertem";
            this.send_libertem.Size = new System.Drawing.Size(50, 22);
            this.send_libertem.TabIndex = 64;
            this.send_libertem.Text = "send";
            this.send_libertem.UseVisualStyleBackColor = true;
            this.send_libertem.Click += new System.EventHandler(this.send_libertem_Click);
            // 
            // chs_blocked
            // 
            this.chs_blocked.AutoSize = true;
            this.chs_blocked.Location = new System.Drawing.Point(330, 172);
            this.chs_blocked.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.chs_blocked.Name = "chs_blocked";
            this.chs_blocked.Size = new System.Drawing.Size(97, 17);
            this.chs_blocked.TabIndex = 65;
            this.chs_blocked.Text = "CH0-6 blocked";
            this.chs_blocked.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            this.chs_blocked.UseVisualStyleBackColor = true;
            this.chs_blocked.CheckedChanged += new System.EventHandler(this.chs_blocked_CheckedChanged);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(856, 484);
            this.Controls.Add(this.chs_blocked);
            this.Controls.Add(this.send_libertem);
            this.Controls.Add(this.libertem_UI);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.label17);
            this.Controls.Add(this.maxNumSlices);
            this.Controls.Add(this.label19);
            this.Controls.Add(this.Download4D);
            this.Controls.Add(this.PreName);
            this.Controls.Add(this.arinaON);
            this.Controls.Add(this.LowResScanTime);
            this.Controls.Add(this.label16);
            this.Controls.Add(this.zerobeam);
            this.Controls.Add(this.alignDif);
            this.Controls.Add(this.magnification_select);
            this.Controls.Add(this.CheckMag);
            this.Controls.Add(this.label15);
            this.Controls.Add(this.BiasP);
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
            this.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.Name = "Form1";
            this.Text = "Form1";
            this.Load += new System.EventHandler(this.Form1_Load);
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
            ((System.ComponentModel.ISupportInitialize)(this.BiasP)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.magnification_select)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.LowResScanTime)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.maxNumSlices)).EndInit();
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
        public System.Windows.Forms.NumericUpDown BiasP;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.Button CheckMag;
        private System.Windows.Forms.NumericUpDown magnification_select;
        private System.Windows.Forms.Button alignDif;
        private System.Windows.Forms.Button zerobeam;
        private System.Windows.Forms.Label label16;
        private System.Windows.Forms.Button arinaON;
        private System.Windows.Forms.TextBox PreName;
        private System.Windows.Forms.Button Download4D;
        private System.Windows.Forms.Label label19;
        public System.Windows.Forms.NumericUpDown LowResScanTime;
        private System.Windows.Forms.NumericUpDown maxNumSlices;
        private System.Windows.Forms.Label label17;
        public System.Windows.Forms.ComboBox libertem_UI;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Button send_libertem;
        private System.Windows.Forms.CheckBox chs_blocked;
    }
}

