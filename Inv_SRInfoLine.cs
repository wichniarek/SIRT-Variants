using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc: Every Singal From S to R
********************************************/

namespace ClassLibrary_TomoGo
{
    public class Inv_SRInfoLine
    {
        public Inv_SRInfoLine()
        {
            name = "Default Name";
            start = new Inv_SR3D();
            end = new Inv_SR3D();
            travelTime = 0f;
        }
        //*********************************************
        private string name;
        public string Name
        {
            get { return name; }
            set { name = value; }
        }
        //*********************************************
        private Inv_SR3D start;
        public Inv_SR3D Start
        {
            get { return start; }
            set { start = value; }
        }
        //---------------------------------------------
        private Inv_SR3D end;
        public Inv_SR3D End
        {
            get { return end; }
            set { end = value; }
        }
        //---------------------------------------------
        private float travelTime;
        public float TravelTime
        {
            get { return travelTime; }
            set { travelTime = value; }
        }
        //---------------------------------------------
        //*********************************************
        public Inv_SRInfoLine copy()
        {
            Inv_SRInfoLine output = new Inv_SRInfoLine();
            output.name = string.Copy(this.name);
            output.start = this.start.copy();
            output.end = this.end.copy();
            output.travelTime = this.travelTime;
            return output;
        }
        //---------------------------------------------
        public string Print(int n)
        {
            //parameter 0: print title and daten
            //parameter 1: only print daten
            string str = null;
            if (n != 0 && n != 1)
            {
                //MessageBox.Show("Function parameter error.");
            }
            else
            {
                str = name + " "
                    + start.Name + "    "
                    + end.Name + " "
                    + travelTime.ToString() + "\r\n";
                if (n == 0)
                {
                    str = "Name" + "    " + "Start" + "  " + "End" + "  " + "TravelTime" + "  " + "\r\n" + str;
                }
                else
                {
                    str = null + str;
                }
            }
            return str;
        }
        //---------------------------------------------
    }
}
