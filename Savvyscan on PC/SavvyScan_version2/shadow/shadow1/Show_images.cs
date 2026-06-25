using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Structure;
using Emgu.CV.Util;
using System;
using System.Drawing;

//using System.Runtime.InteropServices;

namespace shadow1
{
    public class align
    {
        public static bool find_pixel_shift(double angle, ref double shiftx_um, ref double shifty_um)
        {
            //if positive or negative load lastHR_positive.mrc or negative file, and compare with 4 ROIs in the fast scan LR4alignment.mrc, assuming they cover the same area size
            MrcStack myMrcStack = new MrcStack();
            Mat[] slices;
            int shifty_px = 0, shiftx_px = 0;
            int Nrows, Nrows2;
            int Ncols, Ncols2;
            double[] minVal, maxVal;
            Point[] minLoc, maxLoc;

            string FileName = Program.immediate_folder + "newLR4alignment.mrc"; //Low resolution last scan at new tilt angle
            slices = myMrcStack.FOpen(FileName);
            Nrows = slices[0].Rows;
            Ncols = slices[0].Cols;
            Mat LR_mat = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat LR_mat8U = new Mat(Nrows, Ncols, DepthType. Cv8U, 1);
            slices[0].ConvertTo(LR_mat, DepthType.Cv32F, 1, 0);
            CvInvoke.GaussianBlur(LR_mat, LR_mat, new Size(0,0), 6);
            LR_mat.MinMax(out minVal, out maxVal, out minLoc, out maxLoc);
            if (maxVal[0] - minVal[0] > 0) 
                CvInvoke.ConvertScaleAbs(LR_mat, LR_mat8U, 127 / (maxVal[0] - minVal[0]), -minVal[0] * (127 / (maxVal[0] - minVal[0])));
            string FileName2;
            if (angle >= 0)
                FileName2 = Program.immediate_folder + "lastHR_positive.mrc"; //High resolution scan previously collected at near tilt
            else
                FileName2 = Program.immediate_folder + "lastHR_negative.mrc"; //High resolution scan previously collected at near tilt
            slices = myMrcStack.FOpen(FileName2);
            Nrows2 = slices[0].Rows;
            Ncols2 = slices[0].Cols;
            Mat HR_mat = new Mat(Nrows2, Ncols2, DepthType.Cv32F, 1);
            Mat HR2LR_mat = new Mat(Nrows, Ncols, DepthType.Cv32F, 1);
            Mat HR2LR_mat8U = new Mat(Nrows, Ncols, DepthType.Cv8U, 1);
            slices[0].ConvertTo(HR_mat, DepthType.Cv32F, 1, 0);
            Size frame = new Size(Ncols, Nrows);//width,height
            CvInvoke.GaussianBlur(HR_mat, HR_mat, new Size(0, 0), 6);
            CvInvoke.ResizeForFrame(HR_mat, HR2LR_mat, frame, Inter.Linear);
            HR2LR_mat.MinMax(out minVal, out maxVal, out minLoc, out maxLoc);
            if (maxVal[0] - minVal[0]>0)
                CvInvoke.ConvertScaleAbs(HR2LR_mat, HR2LR_mat8U, 127/(maxVal[0]- minVal[0]), -minVal[0]*(127 / (maxVal[0] - minVal[0])));
            Image<Gray, byte> last_image = HR2LR_mat8U.ToImage<Gray, byte>(false);

 
            Mat bottom_roy = new Mat(LR_mat8U, new Rectangle((int)(Ncols / 3), (int)(2 * Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner ,width,height 
            Mat upper_roy = new Mat(LR_mat8U, new Rectangle((int)(Ncols / 3), (int)(0 * Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner:x,y ,width,height 
            Mat right_roy = new Mat(LR_mat8U, new Rectangle((int)(2 * Ncols / 3), (int)(Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner ,width,height 
            Mat left_roy = new Mat(LR_mat8U, new Rectangle((int)(0 * Ncols / 3), (int)(Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner ,width,height 
            Mat br_roy = new Mat(LR_mat8U, new Rectangle((int)(2*Ncols / 3), (int)(2 * Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner ,width,height 
            Mat bl_roy = new Mat(LR_mat8U, new Rectangle((int)(0*Ncols / 3), (int)(2 * Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner:x,y ,width,height 
            Mat tr_roy = new Mat(LR_mat8U, new Rectangle((int)(2 * Ncols / 3), (int)(0*Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner ,width,height 
            Mat tl_roy = new Mat(LR_mat8U, new Rectangle((int)(0 * Ncols / 3), (int)(0*Nrows / 3), (int)(Ncols / 3), (int)(Nrows / 3)));//upper-left corner ,width,height 
            double max = 0;
            Point shift=new Point(0,0);
            Gray avg;
            MCvScalar sdv,sdvref;
            int x0, y0;
            last_image.AvgSdv(out avg, out sdvref);
            for (int n = 1; n <= 8; n++)
            {
                Image<Gray, byte> image_segment;
                switch (n)
                {
                    case 1:
                        image_segment = bottom_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(Ncols / 3); //match result is the new top-left corner location of the segment in the last HR2LR image 
                        y0 = (int)(2 * Nrows / 3);
                        break;
                    case 2:
                        image_segment = upper_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(Ncols / 3);
                        y0 = (int)(0 * Nrows / 3);
                        break;
                    case 3:
                        image_segment = right_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(2 * Ncols / 3);
                        y0 = (int)(Nrows / 3);
                        break;
                    case 4:
                        image_segment = left_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(0 * Ncols / 3);
                        y0 = (int)(Nrows / 3);
                        break;
                    case 5:
                        image_segment = br_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(2* Ncols / 3); //match result is the new top-left corner location of the segment in the last HR2LR image 
                        y0 = (int)(2 * Nrows / 3);
                        break;
                    case 6:
                        image_segment = bl_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(0* Ncols / 3);
                        y0 = (int)(2 * Nrows / 3);
                        break;
                    case 7:
                        image_segment = tr_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(2 * Ncols / 3);
                        y0 = (int)(0*Nrows / 3);
                        break;
                    case 8:
                        image_segment = tl_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(0 * Ncols / 3);
                        y0 = (int)(0 * Nrows / 3);
                        break;
                    default:
                        image_segment = bottom_roy.ToImage<Gray, byte>(false);
                        x0 = (int)(Ncols / 3);
                        y0 = (int)(2 * Nrows / 3);
                        break;
                }
                image_segment.AvgSdv(out avg,out sdv);
                if ((double)sdv.V0 < 0.3 * (double)sdvref.V0) continue; //don't correlate vaccum to vacuum
                var result = last_image.MatchTemplate(image_segment, TemplateMatchingType.CcoeffNormed);//search the segment image (called template) found in argument one in the main image
                result.MinMax(out minVal, out maxVal, out minLoc, out maxLoc);
                if (maxVal[0] > max)
                {
                    max = maxVal[0];
                    shift = maxLoc[0]; 
                    shifty_px = shift.X-x0; //note exchaning x and y and - all based on loading the mrc file to SerialEM and watching the actual shift of the new LR compared to previous stored image 
                    shiftx_px = -(shift.Y-y0);
                }
            }
            shiftx_um = ((double)shiftx_px / Ncols) * (double)Program.LRcella * 0.0001; //convert to microns (cella is standard size of image in angstroms).
            //found that the imageshift is assumed in tilted plane so we don't compensate for cos:   shifty_um = ((double)shifty_px / Nrows) * (double)Program.LRcella * 0.0001/ Program.aspectR; //convert to microns (cella is standard size of image in angstroms).
            //after experiment on 22Apr22, I understood that shift should be in actual pixel size, so it is dy0*cos(theta)=dy0*aspect_ratio
            shifty_um = ((double)shifty_px / Nrows) * (double)Program.LRcella * 0.0001*Program.aspectR; //convert to microns (cella is standard size of image in angstroms).
            if (max > 0.5)
                return true;
            else
                return false;
        }
        private bool Detect_objects(Image<Gray, Byte> Input_Image, Image<Gray, Byte> object_Image)
        {
            //obtined from juanluislm/TeamVis, remove
            Point dftSize = new Point(Input_Image.Width + (object_Image.Width * 2), Input_Image.Height + (object_Image.Height * 2));
            bool Success = false;
            using (Image<Gray, Byte> pad_array = new Image<Gray, Byte>(dftSize.X, dftSize.Y))
            {
                //copy centre
                pad_array.ROI = new Rectangle(object_Image.Width, object_Image.Height, Input_Image.Width, Input_Image.Height);
                CvInvoke.cvCopy(Input_Image.Convert<Gray, Byte>(), pad_array, IntPtr.Zero);
                // CvInvoke.cvMatchTemplate
                //CvInvoke.cvShowImage("pad_array", pad_array);
                pad_array.ROI = (new Rectangle(0, 0, dftSize.X, dftSize.Y));
                using (Image<Gray, float> result_Matrix = pad_array.MatchTemplate(object_Image, TemplateMatchingType.CcoeffNormed))
                {
                    result_Matrix.ROI = new Rectangle(object_Image.Width, object_Image.Height, Input_Image.Width, Input_Image.Height);

                    Point[] MAX_Loc, Min_Loc;
                    double[] min, max;
                    result_Matrix.MinMax(out min, out max, out Min_Loc, out MAX_Loc);

                    using (Image<Gray, double> RG_Image = result_Matrix.Convert<Gray, double>().Copy())
                    {
                        //#TAG WILL NEED TO INCREASE SO THRESHOLD AT LEAST 0.8...used to have 0.7

                        if (max[0] > 0.85)
                        {
                            //Object_Location = MAX_Loc[0];
                            Success = true;
                        }
                    }

                }
            }
            return Success;
        }
    }
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
