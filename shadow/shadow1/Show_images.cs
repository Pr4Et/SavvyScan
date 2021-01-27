using System;
using System.IO;
using Emgu.CV;
using Emgu.Util;
using Emgu.CV.Util;
using Emgu.CV.Structure;
using Emgu.CV.CvEnum;
using System.Drawing;
using System.Windows.Forms;

//using System.Runtime.InteropServices;

namespace shadow1
{
    public class showimage
    {
        public static void Go()
        {
            // path is usually "D:\\ShadowData\\";
            MrcStack myMrcStack = new MrcStack();
            Mat[] slices;
            int Nslices;
            int nslice = 0;
            string win1="", prevwin="";
            int Nrows;
            int Ncols;
            int min, max;
            double minmax = 0;
            double scale;

            for (int ch=0; ch<8;ch++)
            {
                string FileName = Program.immediate_folder + "CH"+ch.ToString()+".mrc"; //dataset, may be tif or mrc
                slices = myMrcStack.FOpen(FileName);
                Nslices = slices.Length;
                win1 = "CH"+ch.ToString();
                Nrows = slices[0].Rows;
                Ncols = slices[0].Cols;
                Mat slice_mat = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
                //Mat slice_mat_win = new Mat(900, 900, DepthType.Cv32F, 1);
                slices[nslice].ConvertTo(slice_mat, DepthType.Cv32F, 1, 0);
                minmax = FThreshold(slice_mat);
                min = (int)minmax;
                max = (int)(100000.0 * (minmax - (double)min));
                scale =255.0/(max-min);
                CvInvoke.ConvertScaleAbs(slice_mat, slice_mat, scale, -scale * min);
                //CvInvoke.Resize(slice_mat, slice_mat_win, new Size(900,900));
                if (ch > 0)
                { CvInvoke.DestroyWindow(prevwin); }

                CvInvoke.NamedWindow(win1, NamedWindowType.FreeRatio); //Create the window using the specific name
                //CvInvoke.LUT(slice_mat,lut, slice_mat) ;
                CvInvoke.Imshow(win1, slice_mat); //Show the image
                CvInvoke.NamedWindow(win1, NamedWindowType.Fullscreen);
                CvInvoke.WaitKey(Program.glancedelay);  //Wait for time/ key pressing event
                prevwin = win1;

            }
            CvInvoke.DestroyWindow(win1);
        }

        static double FThreshold(Mat source)
        {
            double minVal = 0;
            double maxVal = 0;
            System.Drawing.Point locationh = new System.Drawing.Point(0, 0);
            System.Drawing.Point locationl = new System.Drawing.Point(0, 0);
            CvInvoke.MinMaxLoc(source, ref minVal, ref maxVal, ref locationl, ref locationh);
            int size = 4096;
            int[] channels = new int[] { 1 };
            int[] histsize = new int[] { size };
            Mat hist = new Mat(1, size, DepthType.Cv16S, 1);
            Array arr_hist = new int[size, 1];
            //IInputArray mask = null; //null means to ignore
            float[] histrange = new float[] { (float)minVal, (float)maxVal };
            VectorOfMat sources = new VectorOfMat(source);
            sources.Push(source);
            bool accumulate = false;
            CvInvoke.CalcHist(sources, channels, null, hist, histsize, histrange, accumulate);
            arr_hist = hist.GetData();  //from mat to array
            int minInd = 0;
            int maxInd = 0;
            int sum=0;
            for (int ind = 0; ind < size; ind++)
            {
                if (Convert.ToInt32(arr_hist.GetValue(ind, 0)) > 0)
                {
                    sum += Convert.ToInt32(arr_hist.GetValue(ind, 0));
                }
            }
            int sumcomp = 0;
            for (int ind = 0; ind < size; ind++)
            {
                if (Convert.ToInt32(arr_hist.GetValue(ind, 0)) > 0)
                {
                    sumcomp += Convert.ToInt32(arr_hist.GetValue(ind, 0));
                    if (sumcomp/(1.0+sum)>0.002 && minInd==0)
                        minInd = ind;
                    if ((sum-sumcomp) / (1.0 + sum) < 0.002)
                    {
                        maxInd = ind;
                    }
                }
            }
            int min;
            int max;
            min = (int)(minVal + minInd * (maxVal - minVal) / size);
            max = (int)(minVal + maxInd * (maxVal - minVal) / size);
            if (max < min) max = min+1;
            return (double)(min + max / 100000.0);
        }
    }



}
