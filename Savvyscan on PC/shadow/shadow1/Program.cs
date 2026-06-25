using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Text;
using System.IO;

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
        public static string endpoint = Shadow1.Default.serverIP+":5555"; // in lab:"tcp://132.77.57.159:5555", in EM: "tcp://10.0.1.5:5555"
        //public static string endpoint2 = Shadow1.Default.serverIP + ":5556"; 
        public static bool end_request = false;
        public static bool monitor_file=false;
        public static string file_write;
        public static string file_read;
        public static string file_read_defocus= @"C:\Users\stem\Lothar\defocus.txt";
        //public static string file_log = @"D:\SavvyscanData\logfile.txt";
        public static string file_log = @"D:\SavvyscanData\logfile.txt";
        public static string command_from_SerialEM = "";
        static Form1 frm1;
        static ZContext context;
        public static bool busy_talking = false; 
        public static ZSocket requester; //as a client I make requests
        //public static ZSocket responder; //to fetch data as a client I need to respond
        public static ZError error;
        public static int[] data_array= { 0, 0, 0, 0,0 };
        public static string line_defocus;
        public static string immediate_folder="";
        public static int glancedelay = 500;
        public static bool activate_log_file;
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
            // Connect for channel 5555: I talk server replies
            requester = new ZSocket(context, ZSocketType.REQ);
            requester.Connect(endpoint);
            // Connect for channel 5556: server talks I reply
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
            sendrecieve( "outputampmv="+ Shadow1.Default.defaultoutputamp.ToString());
            frm1.OutputScanAmp.Value = Shadow1.Default.defaultoutputamp ;
            sendrecieve( "inputampmv="+ Shadow1.Default.defaultinputamp.ToString());
            frm1.InputScanAmp.Value = Shadow1.Default.defaultinputamp;
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
                            element=word.Split('=');
                            data_array[0] = 1;
                            if (element[0] == "ch1")
                                data_array[1] = int.Parse(element[1]);
                            if (element[0] == "ch2")
                                data_array[2] = int.Parse(element[1]);
                            if (element[0] == "ch3")
                                data_array[3] = int.Parse(element[1]);
                            if (element[0] == "ch4")
                                data_array[4] = int.Parse(element[1]);
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
            busy_talking = false;
            return 0;

        }
        public static string ReplaceInvalidChars(string filename)
        {
            return string.Join("_", filename.Split(Path.GetInvalidFileNameChars()));
        }
    }
}
