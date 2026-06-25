using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace shadow1
{
    public partial class Form1 : Form
    {
        Bitmap DrawArea;
        public Form1()
        {
            InitializeComponent();
            DrawArea = new Bitmap(pictureBox1.Size.Width, pictureBox1.Size.Height);
            pictureBox1.Image = DrawArea;
            Program.file_write =Folder2Monitor.Text + @"\defocus.txt";
            Program.file_read = Folder2Monitor.Text + @"\Command2Shadow.txt";
            //glancefiles_folder.Text = "D:\\SavvyscanData\\";
            glancefiles_folder.Text = "D:\\SavvyscanData\\";
            Program.immediate_folder = glancefiles_folder.Text;

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
            string targetPath = @"d:\SavvyscanData\" + Program.ReplaceInvalidChars(FolderName.Text);
            string text = "save="+ targetPath;
            Program.sendrecieve(text);
        }

        private void OutputScanAmp_ValueChanged(object sender, EventArgs e)
        {
            string text = "Outputampmv=" + OutputScanAmp.Value;
            Program.sendrecieve(text);
            magnification_combo_SelectedIndexChanged(sender, e);
        }

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

        private void timer1_Tick(object sender, EventArgs e)
        {
            if (!Program.busy_talking)
                Program.fetch_data();
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
            if (Program.monitor_file || checkBox1.Checked==true)
            {
                StreamReader sr_defocus = new StreamReader(Program.file_read_defocus);
                Program.line_defocus = sr_defocus.ReadLine();
                Defocus_view.Text = Program.line_defocus;
                sr_defocus.Close();
            }

            if (Program.monitor_file)
            {
                //copy defocus file from Lothar's folder to the Semur folder
                Program.file_write = Folder2Monitor.Text + @"\defocus.txt";
                 StreamWriter sw = new StreamWriter(Program.file_write);
                sw.WriteLine(Program.line_defocus);
                sw.Close();

                //Get commands from the Semur folder related to SerialEM script
                Program.file_read = Folder2Monitor.Text + @"\Command2Shadow.txt";
                StreamReader sr = new StreamReader(Program.file_read);
                string line = sr.ReadLine();
                sr.Close();
                if (!line.Equals(Program.command_from_SerialEM))//if the command is new
                {
                    Program.command_from_SerialEM = line;
                    if (line.Length>=4)
                    {
                        if (line.Equals("tomography"))
                        {
                            startTomo_Click(sender, e);
                            return;
                        }
                        if (line.Substring(0,4).Equals("next"))
                        {
                            NextFrame_Click(sender, e);//consider fetching the tilt angle from after the 4 letters
                            return;
                        }
                        if (line.Equals("tomoend"))
                        {
                            EndTomo_Click(sender, e);
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
            int total=Program.data_array[1]+ Program.data_array[2]+ Program.data_array[3]+ Program.data_array[4]+4*32767;
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
            //g.FillEllipse(mybrush2,x+x0-R,y+y0-R,R+R,R+R);
            g.DrawEllipse(mypen2, x + x0-R, y + y0-R , 2*R, 2*R);// define the rectangualar corner and widths of containing the ellipse
            g.DrawLine(mypen, x0, y0-15, x0, y0+15);
            g.DrawLine(mypen, x0-15, y0, x0+15, y0);
            pictureBox1.Image = DrawArea;
            g.Dispose();
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
            if (TiltAngle.Value==0)
                TiltAngle.Value = -60;
            AspectRatio.Value = Convert.ToDecimal(Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0));
            string text2 = "tomography=0" ;
            Program.sendrecieve(text2);
            int numofslices = 1 + Convert.ToInt16(Math.Abs(Convert.ToDouble(TiltAngle.Value))*2/ Convert.ToDouble(StepAngle.Value));
            string text3 = "numberOfSlices=" + numofslices.ToString();
            Program.sendrecieve(text3);

        }

        private void NextFrame_Click(object sender, EventArgs e)
        {
            tomogramIndex.Value += 1;
            TiltAngle.Value += StepAngle.Value;
            AspectRatio.Value = Convert.ToDecimal(Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0));
            string text2 = "tomography="+ tomogramIndex.Value.ToString();
            Program.sendrecieve(text2);

        }

        private void EndTomo_Click(object sender, EventArgs e)
        {
            TiltAngle.Value = 0;
            tomogramIndex.Value = 0;
            AspectRatio.Value = 1;
            string text2 = "tomography=-1";//sign to end tomography session
            Program.sendrecieve(text2);
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
                checkBox1.Checked = true;
            }
            else
            {
                MonitorFile.BackColor = default(Color);
            }
            System.Threading.Thread.Sleep(200);
        }

        private void magnification_combo_SelectedIndexChanged(object sender, EventArgs e)
        {
            double magnification=9900;
            double.TryParse(magnification_combo.SelectedItem.ToString(), out magnification);
            //based on calibration_17sep20
            int cella = (int)((43628.0 / magnification) * 10.0 * 2048.0 * ((double)OutputScanAmp.Value / 2000.0)); //size of image in angstorms (goes to mrc file, cella.x and cella.y) 
            string text2 = "cella=" + cella.ToString();
            Program.sendrecieve(text2);

        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox1.Checked==true)
            {
                string text2 = "Lothar=1";
                Program.sendrecieve(text2);
            }
            else
            {
                string text2 = "Lothar=0";
                Program.sendrecieve(text2);
            }
        }

        private void TiltAngle_ValueChanged(object sender, EventArgs e)
        {
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
            StartMulti.BackColor = Color.Crimson;
            TiltAngle.Value = 0;
            StepAngle.Value = 0;
            AspectRatio.Value = Convert.ToDecimal(Math.Cos(Convert.ToDouble(TiltAngle.Value) * Math.PI / 180.0));
            string text2 = "tomography=0";
            Program.sendrecieve(text2);
            int numofslices = 300;
            string text3 = "numberOfSlices=" + numofslices.ToString();
            Program.sendrecieve(text3);

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
