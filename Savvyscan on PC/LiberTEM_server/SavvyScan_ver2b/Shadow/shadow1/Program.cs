using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Text;
using System.IO;
using System.Diagnostics;

//using NetMQ;
//using NetMQ.Sockets;
//using Microsoft.SqlServer.Server;
using ZeroMQ;
using System.ServiceModel;

namespace shadow1
{
    public class Program
    {
        //connect to server, fixed IP
        public static string endpoint = Shadow1.Default.serverIP+":5555"; //connect to camera server 
        public static string endpoint2 = Shadow1.Default.serverIP + ":5580"; //connect to python service
        public static bool end_request = false;
        public static bool monitor_file=false;
        public static bool ArinaON = false;
        public static string file_write;
        public static string file_read;
        public static string file_read_defocus= @"C:\Users\stem\Lothar\defocus.txt";
        //public static string file_log = @"D:\SavvyscanData\logfile.txt";
        public static string file_log = @"D:\SavvyscanData\logfile.txt";
        public static string command_from_SerialEM = "";
        //public static string ArinaRequestRun = @"C:\Users\stem\Savvy\run_Arina.bat";
        //public static string ArinaRequestProto = @"C:\Users\stem\Savvy\Proto_run_Arina.bat";
        public static string ArinaFileDownload = @"C:\Users\stem\Savvy\run_Download.bat";
        static Form1 frm1;
        static ZContext context,context2;
        public static bool busy_talking = false;
        public static bool busy_talking2 = false;
        public static ZSocket requester,requester2; //as a client I make requests
        //public static ZSocket responder; //to fetch data as a client I need to respond
        public static ZError error;
        public static int[] data_array= { 0, 0, 0, 0, 0, 0};
        public static double[] data_array2 = { 0, 0, 0, 0, 0 };
        public static string line_defocus;
        public static string immediate_folder="";
        public static int glancedelay = 500;
        public static bool activate_log_file;
        public static int cella;
        public static int LRcella;
        public static double aspectR=1;
        public static bool multislice = false;
        public static bool order_slice_advance = false;
        public static bool order_Arinaoff = false;
        public static int BiasFactorP = 0;
        public static double DiffractionShift_x=0;
        public static double DiffractionShift_y=0;
        public static int caldif_index=0;
        public static double[,] caldif = new double[6, 100]; //0-xshift,1-yshift,2-xdisk,3-ydisk,4-ch5,5-2nd_min_of_ch14)
        public static bool flag_dif_alignment=false;
        public static bool flag_search_scan = false;
        public static double DifDiskx = 0;
        public static double DifDisky = 0;
        public static System.Diagnostics.Process process_Arina = new System.Diagnostics.Process();
        public static System.Diagnostics.ProcessStartInfo startInfo_Arina = new System.Diagnostics.ProcessStartInfo();
        public static int threshold_time=8; //threshold to record


        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);

            frm1 = new Form1();
            frm1.ServerAddress.Text=Shadow1.Default.serverIP;
            frm1.Show();
 
            // Create
            context = new ZContext();
            // Connect for channel 5555, SavvyScan camera server for SerialEM: I talk server replies
            requester = new ZSocket(context, ZSocketType.REQ);
            requester.Connect(endpoint);
            context2 = new ZContext();
            // Connect for channel 5580, local python server (for special functions): I talk server replies
            requester2 = new ZSocket(context2, ZSocketType.REQ);
            requester2.Connect(endpoint2);
            // Connect for channel 5557: server talks I reply
            //responder = new ZSocket(context, ZSocketType.REP);
            //responder.Connect(endpoint2);

            string requestText = "Hello";
            frm1.textBox1.AppendText(String.Format("Sending {0}…\n", requestText));

            // Send
            requester.Send(new ZFrame(requestText));
            // Receive
            using ( ZFrame reply = requester.ReceiveFrame())
            {
                        frm1.textBox1.AppendText(String.Format(" Received: {0}\n",  reply.ReadString()));
            }

            sendrecieve( "scanmode="+ Shadow1.Default.defaultscanmode.ToString());
            frm1.ScanMode_comboBox1.SelectedIndex = Shadow1.Default.defaultscanmode-1;
            sendrecieve( "OutputAmpmv="+ Shadow1.Default.defaultoutputamp.ToString());
            frm1.OutputScanAmp.Value = Shadow1.Default.defaultoutputamp ;
            sendrecieve( "InputAmpmv="+ Shadow1.Default.defaultinputamp.ToString());
            frm1.InputScanAmp.Value = Shadow1.Default.defaultinputamp;
            BiasFactorP = (int)frm1.BiasP.Value;
            string text3 = "BiasOutputP=" + Program.BiasFactorP.ToString();
            sendrecieve(text3);
            text3 = String.Format("LowResTimeS={0}", threshold_time);
            Program.sendrecieve(text3);
            //Process process = System.Diagnostics.Process.Start(@"C:\Users\stem\Lothar\run_script2.bat");

            startInfo_Arina.WindowStyle = System.Diagnostics.ProcessWindowStyle.Normal; //Hidden;
            startInfo_Arina.FileName = @"cmd.exe";
            startInfo_Arina.Arguments = "";
            process_Arina.StartInfo = startInfo_Arina;
            process_Arina.StartInfo.UseShellExecute = false;

            System.Diagnostics.Process process = new System.Diagnostics.Process();
            System.Diagnostics.ProcessStartInfo startInfo = new System.Diagnostics.ProcessStartInfo();
            startInfo.WindowStyle = System.Diagnostics.ProcessWindowStyle.Normal; //Hidden;
            startInfo.FileName = @"cmd.exe";
            startInfo.Arguments = "";
            process.StartInfo = startInfo;
            process.StartInfo.UseShellExecute = false;
            process = System.Diagnostics.Process.Start(@"C:\Users\stem\Savvy\run_script2.bat");
            //process.StartInfo.RedirectStandardInput = true;
            //process.StandardInput.WriteLine(@"cd C:\Users\stem\Lothar");
            //process.StandardInput.WriteLine(@"echo DPC-Shift-auto.py");

            //System.Diagnostics.Process.Start(@"C:\Users\stem\Lothar\run_script3.bat");
            //ProcessStartInfo start = new ProcessStartInfo();
            //start.FileName = @"C:\Windows\System32\cmd.exe";
            //start.WorkingDirectory = @"D:\script";
            //start.Arguments = @"/K C:\Users\stem\Lothar\run_script3.bat";
            //start.UseShellExecute = false;
            //start.RedirectStandardOutput = true;
            //Process process = Process.Start(start);

            /*using (Process process = System.Diagnostics.Process.Start(start))
            {
                using (StreamReader reader = process.StandardOutput)
                {
                    string result = reader.ReadToEnd();
                    Console.Write(result);
                }
            }
            */

            Application.Run(frm1);

            //while (!end_request) { }
            //System.Threading.Thread.Sleep(0);


        }
        public static int sendrecieve(string command)
        {
            busy_talking = true;
            requester.Send(new ZFrame(command));
            using (ZFrame reply = requester.ReceiveFrame())
            {
                if (reply.ReadString() == "OK")
                {
                    busy_talking = false;
                    return 0;
                }
                else
                    return -1;
            }
        }
        public static int fetch_data()
        {
            busy_talking = true;
            requester.Send(new ZFrame("data?"));
            using (ZFrame reply = requester.ReceiveFrame())
            {
                string text = reply.ReadString();
                if (text=="none")
                {
                    data_array[0] = 0;

                }
                else
                {
                    frm1.textBox1.AppendText(text + "\n");
                    string[] words = text.Split(',');
                    string[] element;
                    foreach (var word in words)
                    {
                        if (word.Length>0)
                        {
                            if (word=="Use for alignment")
                            {
                                frm1.align_shift();
                            }
                            else if (word == "Use for @SerialEM")
                            {
                                //do nothing. image is for internal SerialEM use
                            }
                            else
                            { 
                                element=word.Split('=');
                                data_array[0] = 1;
                                if (element[0] == "ch1")
                                    data_array[1] = int.Parse(element[1]);
                                if (element[0] == "ch2")
                                    data_array[2] = int.Parse(element[1]);
                                if (element[0] == "ch3")
                                    data_array[3] = int.Parse(element[1]);
                                if (element[0] == "ch4")
                                {
                                    data_array[4] = int.Parse(element[1]);
                                    if (multislice) order_slice_advance = true;
                                    if (!multislice && ArinaON) order_Arinaoff = true;
                                }
                                if (element[0] == "ch5")
                                    data_array[5] = int.Parse(element[1]) ;
                                if (element[0] == "defocus" && Program.monitor_file)
                                {
                                    float defocus = float.Parse(element[1]);
                                    frm1.Defocus_view.Text = defocus.ToString();
                                    StreamWriter sw = new StreamWriter(file_write);
                                    sw.WriteLine(defocus.ToString());
                                    sw.Close();
                                }
                            }
                        }
                    }

                }
            }
            busy_talking = false;
            return 0;

        }

        public static int sendrecieve_python(string command)
        {
            if (busy_talking2) return -1;
            busy_talking2 = true;
            requester2.Send(new ZFrame(command));
            using (ZFrame reply = requester2.ReceiveFrame())
            {
                if (reply.ReadString() == "OK")
                {
                    busy_talking2 = false;
                    return 0;
                }
                else
                {
                    //busy_talking2 = false;
                    return -1;
                }
            }
        }

        public static int fetch_data_python()
        {
            if (busy_talking2) return -1;
            busy_talking2 = true;
            requester2.Send(new ZFrame("data"));
            using (ZFrame reply = requester2.ReceiveFrame())
            {
                string text = reply.ReadString();
                if (text == "none")
                {
                    data_array2[0] = 0;

                }
                else
                {
                    string[] words = text.Split(',');
                    string[] element;
                    foreach (var word in words)
                    {
                        if (word.Length > 0)
                        {
                                element = word.Split('=');
                                if (element[0] == "x")
                                    data_array2[0] = double.Parse(element[1]);
                                if (element[0] == "y")
                                    data_array2[1] = double.Parse(element[1]);
                         }
                    }

                }
            }
            busy_talking2 = false;
            return 0;

        }

 
        public static string ReplaceInvalidChars(string filename)
        {
            return string.Join("_", filename.Split(Path.GetInvalidFileNameChars()));
        }
    }
}
