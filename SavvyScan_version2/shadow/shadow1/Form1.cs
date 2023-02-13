using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.IO.Ports;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Diagnostics;


namespace shadow1
{
    public partial class Form1 : Form
    {
        Bitmap DrawArea;
        string TextAfter = "";
        string NextTextAfter = "0";
        double[] px=new double[3];
        double[] py=new double[3];
        //internal readonly object ThresholdTime;

        //SerialPort mySerialPort = new SerialPort("COM5", 38400); //Connection to ArinaGate box for trigger to Arina detector recording
        //string ComPort_message="";

        public Form1()
        {
            InitializeComponent();
            DrawArea = new Bitmap(pictureBox1.Size.Width, pictureBox1.Size.Height);
            pictureBox1.Image = DrawArea;
            Program.file_write =Folder2Monitor.Text + @"\defocus.txt";
            Program.file_read = Folder2Monitor.Text + @"\Command2Shadow.txt";
            glancefiles_folder.Text = "A:\\SavvyscanData\\";
            Program.immediate_folder = glancefiles_folder.Text;
            //mySerialPort.DataBits = 8;
            //mySerialPort.StopBits = StopBits.One;
            //mySerialPort.Parity = 0;
            //mySerialPort.Open();
            //mySerialPort.DataReceived += new SerialDataReceivedEventHandler(SerialPortDataReceived);

        }

        private void label3_Click(object sender, EventArgs e)
        {

        }

        private void END_Click(object sender, EventArgs e)
        {
            Program.end_request = true;
            timer1.Stop();
            Program.requester.Close();
            this.Close();
        }

        private void ScanMode_comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            string text = "scanmode=" + (ScanMode_comboBox1.SelectedIndex + 1).ToString();
            Program.sendrecieve(text);
        }

        private void SavetoFolder_Click(object sender, EventArgs e)
        {
            if (tomogramIndex.Value>0) //to close tomogram stack mrc file before save
            {
                EndTomo_Click(sender, e);
                Task.Delay(500);
            }
            
            string targetPath = @"d:\SavvyscanData\" + PreName.Text+"_"+Program.ReplaceInvalidChars(FolderName.Text);
            string text = "save="+ targetPath;
            Program.sendrecieve(text);
        }

        private void OutputScanAmp_ValueChanged(object sender, EventArgs e)
        {
            string text = "Outputampmv=" + OutputScanAmp.Value;
            Program.sendrecieve(text);
            set_cella();
        }
        //private void biasOutput1_ValueChanged(object sender, EventArgs e)
        //{
        //    string text = "BiasOutputmv=" + biasOutput1.Value;
        //    Program.sendrecieve(text);
        //    magnification_combo_SelectedIndexChanged(sender, e);
        //}

        private void InputScanAmp_ValueChanged(object sender, EventArgs e)
        {
            string text = "Inputampmv=" + InputScanAmp.Value;
            Program.sendrecieve(text);
        }

        private void AspectRatio_ValueChanged(object sender, EventArgs e)
        {
            string text = "aspectratio=" + AspectRatio.Value;
            Program.sendrecieve(text);
        }

        private void ScanArgument_ValueChanged(object sender, EventArgs e)
        {
            string text = "scanargument=" + ScanArgument.Value;
            Program.sendrecieve(text);
        }
        //private void SerialPortDataReceived(object sender, SerialDataReceivedEventArgs e)
       // {
       //     ComPort_message = mySerialPort.ReadLine();
       // }

        private void timer1_Tick(object sender, EventArgs e)
        {
            if (!Program.busy_talking)
            {
                Program.fetch_data();
            }
            if (Program.data_array[0] == 1)
            {
                plot_spotter();
                Program.data_array[0] = 0;
                /*if (Program.activate_log_file)
                {
                    StreamWriter swlog = File.AppendText(Program.file_log);
                    swlog.WriteLine(TUIshiftx.Value+" , " + TUIshifty.Value + " , "+ Program.data_array[1].ToString()+" , "+ Program.data_array[2].ToString() + " , " + Program.data_array[3].ToString() + " , " + Program.data_array[4].ToString());
                    swlog.Close();
                }*/
            }
            if (!Program.busy_talking)
            {
                if (Program.order_slice_advance)
                {
                    Program.order_slice_advance = false;
                    NextFrame_Auto();
                }
                if (Program.order_Arinaoff)
                {
                    Program.order_Arinaoff = false;
                    if (Program.ArinaON) arinaON_Click_1(sender, e);//click to remove Arina armining for single scan
                }
                if (Program.flag_search_scan)
                {
                    Program.flag_search_scan = false;
                    Program.sendrecieve_python("sem.Search()");
                }
            }
            if (checkBox1.Checked==true)
            {
                StreamReader sr_defocus = new StreamReader(Program.file_read_defocus);
                Program.line_defocus = sr_defocus.ReadLine();
                Defocus_view.Text = Program.line_defocus;
                sr_defocus.Close();
                if (Program.monitor_file)
                {
                    //copy defocus file from Lothar's folder to the shared folder
                    Program.file_write = Folder2Monitor.Text + @"\defocus.txt";
                    StreamWriter sw = new StreamWriter(Program.file_write);
                    sw.WriteLine(Program.line_defocus);
                    sw.Close();
                }
            }

            if (Program.monitor_file)
            {
 
                //Get commands from the network folder related to SerialEM script
                Program.file_read = Folder2Monitor.Text + @"\Command2Shadow.txt";
                StreamReader sr = new StreamReader(Program.file_read);
                string line = sr.ReadLine();
                sr.Close();
                if (line == null) line = " ";
                if (!line.Equals(Program.command_from_SerialEM))//if the command is new
                {
                    Program.command_from_SerialEM = line;
                    if (line.Length>=4)
                    {
                        if (line.Substring(0,Math.Min(10,line.Length)).Equals("tomography"))
                        {
                            if (line.Length>11)
                            {
                                TextAfter = line.Substring(11, line.Length - 11);
                                FolderName.Text = TextAfter;
                                StartMulti_Click_1(sender, e);
                            }
                            return;
                        }
                        if (line.Substring(0, Math.Min(9, line.Length)).Equals("startTomo"))
                        {
                            if (line.Length > 10)
                            {
                                TextAfter = line.Substring(10, line.Length - 10);
                                FolderName.Text = TextAfter;
                                startTomo_Auto();
                            }
                            return;
                        }
                        if (line.Substring(0,4).Equals("next") && !Program.multislice)
                        {
                            if (line.Length > 5)
                            {
                                NextTextAfter = line.Substring(5, line.Length - 5);
                            }
                            else
                            {
                                NextTextAfter = "";
                            }
                            NextFrame_Auto();//for fetching the tilt angle from after the 4 letters
                            return;
                        }
                        if (line.Equals("tomoend"))
                        {
                            EndTomo_Click(sender, e);
                            SavetoFolder_Click(sender, e); //Save automaticaly to the folder written in GUI
                            return;
                        }
                     }
                }
            }
        }

        public void plot_spotter()
        {
            Graphics g;
            g = Graphics.FromImage(DrawArea);
            Pen mypen = new Pen(Color.Black);
            Pen mypen2 = new Pen(Color.Red);
            //Brush mybrush2= new Brush(Color.Red);
            g.Clear(Color.Aqua);
            //int total = Program.data_array[1] + Program.data_array[2] + Program.data_array[3] + Program.data_array[4];// +4*32767;
            float x0 = DrawArea.Width / 2;
            float y0 = DrawArea.Height / 2;
            float x = 0;
            float y = 0;
            float R = 40;
            //calculate from signal interpreted as areas of circle sections 
            float sumup = Program.data_array[1] + Program.data_array[4];
            float sumdown = Program.data_array[2] + Program.data_array[3];
            float sumright = Program.data_array[1] + Program.data_array[2];
            float sumleft = Program.data_array[4] + Program.data_array[3];
            float theta1 = solve_theta(Math.PI*(1+Math.Abs(sumdown- sumup)/(Math.Abs(sumdown)+ Math.Abs(sumup) + 100)));
            float theta2 = solve_theta(Math.PI *(1+ Math.Abs(sumright - sumleft) / (Math.Abs(sumright) + Math.Abs(sumleft) + 100)));
            x = R * (float)Math.Sqrt(1.0 - Math.Pow(Math.Sin(theta2/2), 2));
            y = R * (float)Math.Sqrt(1.0 - Math.Pow(Math.Sin(theta1/2), 2));
            //mind the signs
            if (sumup<sumdown) x=-x;
            if (sumright < sumleft) y = -y;
            Program.DifDiskx = x;
            Program.DifDiskx = y;
            //g.FillEllipse(mybrush2,x+x0-R,y+y0-R,R+R,R+R);
            g.DrawEllipse(mypen2, x + x0-R, y + y0-R , 2*R, 2*R);// define the rectangualar corner and widths of containing the ellipse
            g.DrawLine(mypen, x0, y0-15, x0, y0+15);
            g.DrawLine(mypen, x0-15, y0, x0+15, y0);
            pictureBox1.Image = DrawArea;
            g.Dispose();
            if (Program.flag_dif_alignment) //in mode=6 (dif shift alignement)
            {
                //Fill line in caldif table
                //0-xshift,1-yshift,2-xdisk,3-ydisk,4-ch5,5-2nd_min_of_ch14)
                double cursor_throld = 15;
                Program.caldif[0, Program.caldif_index] = Program.DiffractionShift_x;
                Program.caldif[1, Program.caldif_index] = Program.DiffractionShift_y;
                Program.caldif[2, Program.caldif_index] = x; //xdisk position
                Program.caldif[3, Program.caldif_index] = y; //ydisk position
                Program.caldif[4, Program.caldif_index] = Program.data_array[5]; //ch5
                double[] data2sort = new double[4];
                for (int ind = 0; ind < 4; ind++) { data2sort[ind] = Program.data_array[ind+1]; };
                Array.Sort(data2sort);
                Program.caldif[5, Program.caldif_index] = data2sort[0]; //The least value of ch1 to 4
                if (Program.caldif_index < 100 - 1) Program.caldif_index++;
                else { textBox1.AppendText("Too many points\n"); return; }
                //Calibrate parameters and find the best diffraction shift to minimize CH5 signal and to place disk position at 0,0
                int count_vect_max = 0;//determine size of vectors for regression later
                for (int ind = 0; ind < Program.caldif_index; ind++)
                {
                    if (Math.Abs(Program.caldif[2, ind]) < cursor_throld && Math.Abs(Program.caldif[3, ind]) < cursor_throld) //if ch1-4 larger than ch5 and cursor near center Program.caldif[5,ind]> Program.caldif[4, Program.caldif_index] && 
                    {
                        count_vect_max++;
                    }
                }
                double[][] feature = new double[count_vect_max][];
                double[][] feature_ch5 = new double[count_vect_max][];
                double[] labelx = new double[count_vect_max];
                double[] labely = new double[count_vect_max];
                double[] label_ch5 = new double[Program.caldif_index];
                int count_vect = 0;
                double minval = 100000;
                //int min_ind = 0;
                double maxval = 10;
                int best_ind = 0;
                bool succeeded = false;
                bool last_is_in = false;

                bool already_found = false;
                for (int ind=0; ind< Program.caldif_index; ind++)
                {
                    if (Math.Abs(Program.caldif[2, ind])< cursor_throld && Math.Abs(Program.caldif[3, ind]) < cursor_throld) //if ch1-4 larger than ch5 and x,y cursor near center
                    {
                        feature[count_vect] = new[] { Program.caldif[0, ind], Program.caldif[1, ind] };
                        labelx[count_vect] = Program.caldif[2, ind];
                        labely[count_vect] = Program.caldif[3, ind];
                        count_vect++;
                        last_is_in = true;
                    }
                    else
                    {
                        last_is_in = false;
                    }
                    //feature_ch5[ind] = new[] { Program.caldif[0, ind], Program.caldif[1, ind] };
                    //label_ch5[ind] = Program.caldif[4, ind];
                    if (Program.caldif[4, ind]< minval)
                    {
                        minval = Program.caldif[4, ind];
                        best_ind = ind;
                    }
                    if (Program.caldif[5, ind] > maxval && Program.caldif[4, ind]<minval+10 )
                    {
                        maxval = Program.caldif[5, ind];
                        //best_ind = ind;
                        //textBox1.AppendText(String.Format("maxval{0}", maxval));
                    }
                    
                }
                if (Math.Sqrt(x * x + y * y) <= 2 && Program.caldif_index > 3) already_found = true;

                if ((count_vect>=8 || already_found) && last_is_in)
                {
                    px = MathNet.Numerics.LinearRegression.MultipleRegression.QR(feature, labelx, intercept: true); //array double[]
                    py = MathNet.Numerics.LinearRegression.MultipleRegression.QR(feature, labely, intercept: true);
                    var A = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.DenseOfArray(new double[,] {{ px[1], px[2] },{py[1], py[2] }});
                    var b = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(new double[] { -px[0], -py[0] });
                    try
                    {
                        var res = A.Solve(b);
                        Program.DiffractionShift_x = res[0];
                        Program.DiffractionShift_y = res[1];
                        succeeded = true;
                    }
                    catch
                    {
                        Program.DiffractionShift_x = Program.caldif[0, best_ind];
                        Program.DiffractionShift_y = Program.caldif[1, best_ind];
                    }
                }
                else
                {
                    var rand = new Random();
                    Program.DiffractionShift_x = Program.caldif[0, best_ind] + 0.15* rand.NextDouble();
                    Program.DiffractionShift_y = Program.caldif[1, best_ind] + 0.15* rand.NextDouble();
                }
                textBox1.AppendText(String.Format("SetDiffractionShift({0:0.######},{1:0.######})\n", Program.DiffractionShift_x, Program.DiffractionShift_y));
                Program.sendrecieve_python(String.Format("sem.SetDiffractionShift({0:0.######},{1:0.######})", Program.DiffractionShift_x, Program.DiffractionShift_y));
                Task.Delay(500);
                Program.sendrecieve_python("(x,y)=sem.ReportDiffractionShift()");
                Task.Delay(500);
                int feedbk = Program.fetch_data_python();
                if (feedbk != -1)
                {
                    Program.DiffractionShift_x = Program.data_array2[0];
                    Program.DiffractionShift_y = Program.data_array2[1];
                }
             
                if (succeeded)
                {
                    Program.flag_dif_alignment = false;
                    ScanMode_comboBox1.SelectedIndex = 0;
                    Shadow1.Default.p0x = px[0];
                    Shadow1.Default.p1x = px[1];
                    Shadow1.Default.p2x = px[2];
                    Shadow1.Default.p0y = py[0];
                    Shadow1.Default.p1y = py[1];
                    Shadow1.Default.p2y = py[2];
                    textBox1.AppendText(String.Format("px=({0:0.######},{1:0.######},{1:0.######})\n", px[0],px[1],px[2]));
                    textBox1.AppendText(String.Format("py=({0:0.######},{1:0.######},{1:0.######})\n", py[0], py[1], py[2]));
                    alignDif.BackColor = Color.Beige;
                    Program.sendrecieve_python("single"); //end continous mode in python program
                    Task.Delay(3000);
                }
                Task.Delay(500);
                Program.flag_search_scan = true; //request search scan when possible (on clock interrupt)
            }
        }
        private float solve_theta(double C)
        {
            //Solve theta-sin(theta)=C by Newton's method
            double Cnorm = C;
            if (Cnorm > 2*Math.PI) Cnorm = 2*Math.PI;
            if (Cnorm < -2*Math.PI) Cnorm = -2*Math.PI;
            double theta=Cnorm;
            for (int n=0; n<10; n++)
            {
                if (Math.Abs(-1+Math.Cos(theta))>0.01)
                    theta = theta - (Cnorm - theta+Math.Sin(theta)) / (-1 + Math.Cos(theta));
            }
            return (float)theta;

        }

        private void startTomo_Click(object sender, EventArgs e)
        {
            startTomo.BackColor = Color.Crimson;
            Program.multislice = false;
            string text1 = "SerialEMAlign=0";//mark that we are doing the image alignment and not SerialEM
            Program.sendrecieve(text1);
            set_cella();
            if (TiltAngle.Value==0 && StepAngle.Value>0)
                TiltAngle.Value = -60;
            Program.aspectR = Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0);
            AspectRatio.Value = Convert.ToDecimal(Program.aspectR); //will be sensed and send command aspectratio=
            tomogramIndex.Value = 0;
            string text2 = "tomography=0" ;
            Program.sendrecieve(text2);
            int numofslices;
            if (StepAngle.Value>0)
            { numofslices = 1 + Convert.ToInt16(60 * 2 / Convert.ToDouble(StepAngle.Value)); }
            else
            { numofslices = 16; }
            string text3 = "numberOfSlices=" + numofslices.ToString();
            Program.sendrecieve(text3);

        }

        private void NextFrame_Click(object sender, EventArgs e)
        {
            tomogramIndex.Value += 1;
            if (!Program.multislice)
            {
                TiltAngle.Value += StepAngle.Value;
                Program.aspectR = Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0);
                AspectRatio.Value = Convert.ToDecimal(Program.aspectR); //will be sensed and send command aspectratio=
            }
            string text2 = "tomography="+ tomogramIndex.Value.ToString();
            Program.sendrecieve(text2);

        }
        private void startTomo_Auto()
        {
            startTomo.BackColor = Color.Crimson;
            Program.multislice = false;
            string text1 = "SerialEMAlign=0";
            Program.sendrecieve(text1);
            set_cella();
            tomogramIndex.Value = -1;
            string text2 = "tomography=-1"; //-1 so on command next that will send the angle the index will be 0 that is the first slice in tomography
            Program.sendrecieve(text2);
            int numofslices = 121;
            string text3 = "numberOfSlices=" + numofslices.ToString();
            Program.sendrecieve(text3);
            //string text4 = "foldername=" + TextAfter;
            //Program.sendrecieve(text4);
        }

        public void NextFrame_Auto()
        {
            tomogramIndex.Value += 1;
            string text2 = "tomography=" + tomogramIndex.Value.ToString();
            Program.sendrecieve(text2);
            if (!Program.multislice)
            {
                decimal angle = 0;
                decimal.TryParse(NextTextAfter, out angle);
                TiltAngle.Value = angle;//tiltangle change will be sensed, see below (command will be sent to savvyscan)
            }

        }
        public void align_shift()
        {
            double shiftx_um = 0, shifty_um = 0;
            align.find_pixel_shift((double)TiltAngle.Value, ref shiftx_um, ref shifty_um);
            //send alignmet command to SerialEM
            Program.file_write = Folder2Monitor.Text + @"\alignment.txt";
            StreamWriter sw = new StreamWriter(Program.file_write);
            sw.WriteLine((0.0001*Math.Round(10000*shiftx_um)).ToString()+" "+ (0.0001*Math.Round(10000*shifty_um)).ToString());
            sw.Close();

        }
        private void EndTomo_Click(object sender, EventArgs e)
        {
            TiltAngle.Value = 0;
            tomogramIndex.Value = 0;
            Program.aspectR = 1.0;
            AspectRatio.Value = Convert.ToDecimal(1.0);
            string text2 = "tomography=-2";//sign to end tomography session (-2 since must be below -1)
            Program.sendrecieve(text2);
            if (Program.multislice)
            {
                string text1 = "SerialEMAlign=1";
                Program.sendrecieve(text1);
                Program.multislice = false;

                arinaON.BackColor = Color.Beige;
                Program.ArinaON = false;
                text1 = "ArinaON=0";
                Program.sendrecieve(text1);
            }
            startTomo.BackColor = default(Color);
            StartMulti.BackColor = default(Color);
        }

        private void chosenCH_ValueChanged(object sender, EventArgs e)
        {
            string text2 = "ch=" + chosenCH.Value.ToString();
            Program.sendrecieve(text2);

        }

        private void button1_Click(object sender, EventArgs e)
        {
            Program.monitor_file = !Program.monitor_file;
            if (Program.monitor_file)
            {
                MonitorFile.BackColor = Color.Red;
                //checkBox1.Checked = true; //Enable reading defocus from 4-quadrants
            }
            else
            {
                MonitorFile.BackColor = default(Color);
            }
            System.Threading.Thread.Sleep(200);
        }

         private void set_cella()
        {
            double magnification =(double)magnification_select.Value;
            double LRmagnification = 9900;
            //double.TryParse(magnification_combo.SelectedItem.ToString(), out magnification);

            //based on calibration_17sep20
            Program.cella = (int)((43628.0 / magnification) * 10.0 * 2048.0 * ((double)OutputScanAmp.Value / 2000.0)); //size of image in angstorms (goes to mrc file, cella.x and cella.y) 
            Program.LRcella = (int)((43628.0 / LRmagnification) * 10.0 * 2048.0 * ((double)OutputScanAmp.Value / 2000.0)); //size of image in angstorms (goes to mrc file, cella.x and cella.y) 
            string text2 = "cella=" + Program.cella.ToString();
            Program.sendrecieve(text2);
            text2 = "LRcella=" + Program.LRcella.ToString();
            Program.sendrecieve(text2);
        }
        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox1.Checked==true)
            {
                string text2 = "Lothar=1";
                Program.sendrecieve(text2);
                Process cmd = new Process();
                cmd.StartInfo.FileName = @"C:\Windows\System32\cmd.exe";
                cmd.StartInfo.RedirectStandardInput = true;
                cmd.StartInfo.RedirectStandardOutput = true;
                cmd.StartInfo.CreateNoWindow = true;
                cmd.StartInfo.UseShellExecute = false;
                cmd.Start();
                cmd.StandardInput.WriteLine(@"C:\Users\stem\Lothar\run_script.bat");
                cmd.StandardInput.Flush();
                cmd.StandardInput.Close();
                cmd.WaitForExit();
                //Console.WriteLine(cmd.StandardOutput.ReadToEnd());

            }
            else
            {
                string text2 = "Lothar=0";
                Program.sendrecieve(text2);
            }
        }

        private void TiltAngle_ValueChanged(object sender, EventArgs e)
        {
            Program.aspectR = Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0);
            if (Program.multislice) Program.aspectR=1.0;
            AspectRatio.Value = Convert.ToDecimal(Program.aspectR); //will be sensed and send command aspectratio=
            string text4 = "tiltangle=" + TiltAngle.Value.ToString();
            Program.sendrecieve(text4);

        }

        private void saveMAT_CheckedChanged(object sender, EventArgs e)
        {
            if (saveMAT.Checked == true)
            {
                string text2 = "saveMAT=1";
                Program.sendrecieve(text2);
            }
            else
            {
                string text2 = "saveMAT=0";
                Program.sendrecieve(text2);
            }

        }

        private void Save_MATKEY_CheckedChanged(object sender, EventArgs e)
        {
            if (Save_MATKEY.Checked == true)
            {
                if (saveMAT.Checked == false)
                {
                    saveMAT_CheckedChanged(sender, e);
                }
                string text2 = "saveKEY=1";
                Program.sendrecieve(text2);
            }
            else
            {
                string text2 = "saveKEY=0";
                Program.sendrecieve(text2);
            }

        }

        private void GlanceFiles_Click(object sender, EventArgs e)
        {
            showimage.Go();
        }

        private void glancefiles_folder_TextChanged(object sender, EventArgs e)
        {
            Program.immediate_folder = glancefiles_folder.Text;
        }

        private void glance_delay_ValueChanged(object sender, EventArgs e)
        {
            Program.glancedelay = (int)glance_delay.Value;
        }

        private void StartMulti_Click_1(object sender, EventArgs e)
        {
            string text1;
            /*if (Program.ArinaON)
            {
                MessageBox.Show("You must reset Dectris server and arm Arina detector again to take effect");
                arinaON.BackColor = Color.Beige;
                Program.ArinaON = false;
                text1 = "ArinaON=0";
                Program.sendrecieve(text1);
            }*/
            StartMulti.BackColor = Color.Crimson;
            Program.multislice = true;
            text1 = "SerialEMAlign=1";
            Program.sendrecieve(text1);
            set_cella();
            TiltAngle.Value = 0;
            StepAngle.Value = 0;
            Program.aspectR = Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0);
            AspectRatio.Value = Convert.ToDecimal(Program.aspectR); //will be sensed and send command aspectratio=
            tomogramIndex.Value = 0;
            string text2 = "tomography=0";
            Program.sendrecieve(text2);
            int numofslices = (int)maxNumSlices.Value;
            string text3 = "numberOfSlices=" + numofslices.ToString();
            Program.sendrecieve(text3);
            //num_pulses_arm.Value= 1024 * 1024 * 121;
        }

        private void BiasP_ValueChanged(object sender, EventArgs e)
        {
            Program.BiasFactorP = (int)BiasP.Value;
            string text3 = "BiasOutputP=" + Program.BiasFactorP.ToString();
            Program.sendrecieve(text3);
        }

        private void CheckMag_Click(object sender, EventArgs e)
        {
            int feedbk = 0;
            Program.flag_dif_alignment=false;
            Program.caldif_index = 0; //clear caldif table
            Program.flag_search_scan = false;
            alignDif.BackColor = Color.Beige;
            //check magnification
            Task.Delay(500);
            Program.sendrecieve_python("(x,y)=sem.ReportMag()");
            Task.Delay(2000);
            feedbk=Program.fetch_data_python();
            //textBox1.AppendText(String.Format(" Mag= {0}\n", Program.data_array2[0]));
            magnification_select.Value = (decimal)Program.data_array2[0];
            //check Diff Shift
            Program.sendrecieve_python("(x,y)=sem.ReportDiffractionShift()");
            Task.Delay(500);
            feedbk =Program.fetch_data_python();
            if (feedbk == -1) return;
            textBox1.AppendText(String.Format(" DiffractionShift= {0},{1}\n", Program.data_array2[0], Program.data_array2[1]));
            Program.DiffractionShift_x = Program.data_array2[0];
            Program.DiffractionShift_y = Program.data_array2[1];
            Program.sendrecieve_python("single"); //end continous mode in python program
            Task.Delay(500);
        }

        private void magnification_select_ValueChanged(object sender, EventArgs e)
        {
            set_cella();
        }

        private void alignDif_Click(object sender, EventArgs e)
        {
            ScanMode_comboBox1.SelectedIndex = 5;
            ScanMode_comboBox1_SelectedIndexChanged(sender, e);
            alignDif.BackColor = Color.Red;
            Task.Delay(1000);
            Program.sendrecieve_python("continous"); //mode in python program ended by "single"
            Program.flag_dif_alignment = true;
            Task.Delay(1000);
            Program.flag_search_scan = true;
        }

        private void zerobeam_Click(object sender, EventArgs e)
        {
            zerobeam.BackColor= Color.Red;
            string text2 = "CalibrateBias=1";
            Program.sendrecieve(text2);
            Task.Delay(500);
            Program.flag_search_scan = true; //request search scan when possible (on clock interrupt)
            zerobeam.BackColor = Color.Beige;
        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {
            if (Program.DifDiskx == 0 && Program.DifDisky == 0) return;
            px[1] = Shadow1.Default.p1x;
            px[2]= Shadow1.Default.p2x;
            py[1] = Shadow1.Default.p1y;
            py[2] = Shadow1.Default.p2y;

            var A = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.DenseOfArray(new double[,] { { px[1], px[2] }, { py[1], py[2] } });
            //Using only one measurement of dif.shift command and x,y diffraction shift observed in order to find updated bias value, then extract new command to move dif to center.
            px[0] = -(px[1] * Program.DiffractionShift_x + px[2] * Program.DiffractionShift_y - Program.DifDiskx);
            py[0] = -(py[1] * Program.DiffractionShift_x + py[2] * Program.DiffractionShift_y - Program.DifDisky);
            var b = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(new double[] { - px[0], -py[0] });
            try
            {
                var res = A.Solve(b);
                Program.DiffractionShift_x = res[0];
                Program.DiffractionShift_y = res[1];
                textBox1.AppendText(String.Format("SetDiffractionShift({0:0.######},{1:0.######})\n", Program.DiffractionShift_x, Program.DiffractionShift_y));
                Program.sendrecieve_python(String.Format("sem.SetDiffractionShift({0:0.######},{1:0.######})", Program.DiffractionShift_x, Program.DiffractionShift_y));
                Task.Delay(500);
                Program.sendrecieve_python("(x,y)=sem.ReportDiffractionShift()");
                Task.Delay(500);
                int feedbk = Program.fetch_data_python();
                if (feedbk != -1)
                {
                    Program.DiffractionShift_x = Program.data_array2[0];
                    Program.DiffractionShift_y = Program.data_array2[1];
                }
                Program.flag_search_scan = true; //request search scan when possible (on clock interrupt)
            }
            catch
            {
                textBox1.AppendText("|Failed, Try to calibrate diff disk");
            }



        }

        private void numericUpDown1_ValueChanged(object sender, EventArgs e)
        {
            Program.threshold_time = (int)LowResScanTime.Value;
            string text2 = String.Format("LowResTimeS={0}", Program.threshold_time);
            Program.sendrecieve(text2);

        }

        private void arinaON_Click_1(object sender, EventArgs e)
        {
            string text1;
            if (Program.ArinaON)
            { 
                arinaON.BackColor = Color.Beige;
                Program.ArinaON = false;
                text1 = "ArinaON=0";
                Program.sendrecieve(text1);
            }
            else
            {
                arinaON.BackColor = Color.Crimson;
                Program.ArinaON = true;
                text1 = "ArinaON=1";
                Program.sendrecieve(text1);
                string targetPath = PreName.Text +"_"+ Program.ReplaceInvalidChars(FolderName.Text);
                text1 = "ArinaFile=" + targetPath;
                Program.sendrecieve(text1);
            }

            //-n 1048576 -e 200 -o Scan_0 -d c:/ArinaData -x
            //-o PreName.Text+Program.ReplaceInvalidChars(FolderName.Text);
            /*int pulsecount = (int)num_pulses_arm.Value; //expect pulses for one scan normally, 1024X1024, or for entire tilt series in tomography 
            StreamReader sr_proto = new StreamReader(Program.ArinaRequestProto);
            StreamWriter sw_proto = new StreamWriter(Program.ArinaRequestRun);
            string rx = sr_proto.ReadLine();
            sw_proto.WriteLine(rx);
            rx = sr_proto.ReadLine();
            rx=rx.Substring(0, rx.IndexOf(" -n "));
            rx = rx + " -n "+pulsecount.ToString()+" -e 200 -o " + PreName.Text + Program.ReplaceInvalidChars(FolderName.Text)+" -x";//-x is for external trigger! + " -d c:/ArinaData -x";
            sw_proto.WriteLine(rx);
            sr_proto.Close();
            sw_proto.Close();
            Program.process_Arina = System.Diagnostics.Process.Start(Program.ArinaRequestRun);//@"C:\Users\stem\Savvy\run_Arina.bat");
           */

        }

        private void Form1_Load(object sender, EventArgs e)
        {
            PreName.Text = DateTime.Now.ToString("yyyyMMMdd");
        }

        private void Download4D_Click(object sender, EventArgs e)
        {
            Program.process_Arina = System.Diagnostics.Process.Start(Program.ArinaFileDownload);
        }

        private void tomogramIndex_ValueChanged(object sender, EventArgs e)
        {
            string text2 = "tomography=" + tomogramIndex.Value.ToString();
            Program.sendrecieve(text2);

        }




        /*private void LogPosition_Click(object sender, EventArgs e)
        {
            if (Program.activate_log_file==false)
            {
                Program.activate_log_file = true;
                LogPosition.BackColor = Color.Red;
            }
            else
            {
                Program.activate_log_file = false;
                LogPosition.BackColor = default(Color);
            }
        }*/
    }
}
