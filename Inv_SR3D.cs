using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc: The Nodes Class for Inversion
********************************************/

namespace ClassLibrary_TomoGo
{
    public class Inv_SR3D
    {
        public Inv_SR3D()
        {
            name = "Default Name";
            sorR = SorREnum.NotSR; //three kinds of nodes: R, S and others
            coor[0] = float.MaxValue; coor[1] = float.MaxValue; coor[2] = float.MaxValue;//the coordinate of the node
        }
        //-------------------------------------------
        private string name;
        public string Name
        {
            get { return name; }
            set { name = value; }
        }
        //-------------------------------------------
        public enum SorREnum
        {
            Receiver, Sender, NotSR
        }
        private SorREnum sorR;
        public SorREnum SorR
        {
            get { return sorR; }
            set { sorR = value; }
        }
        //-------------------------------------------
        private float[] coor = new float[3];
        public float[] Coor
        {
            get { return coor; }
            set { coor = value; }
        }
        //-------------------------------------------
        private int no = 9999;
        public int No
        {
            get { return no; }
            set { no = value; }
        }
        //-------------------------------------------
        public Inv_SR3D copy()
        {
            Inv_SR3D output = new Inv_SR3D();
            output.name = string.Copy(this.name);
            output.no = this.no;
            output.coor[0] = this.coor[0];
            output.coor[1] = this.coor[1];
            output.coor[2] = this.coor[2];
            return output;
        }
        //-------------------------------------------
        public string print(int n)
        {
            //parameter 0: print title and daten
            //parameter 1: only print daten
            string str = null;
            if (n != 0 && n != 1)
            {
                // MessageBox.Show("Function parameter error.");
            }
            else
            {
                str = name + " "
                    + sorR.ToString() + "    "
                    + coor[0].ToString() + " " + coor[1].ToString() + " " + coor[2].ToString() + " "
                    + no.ToString() + "\r\n";
                if (n == 0)
                {
                    str = "Name" + "    " + "Type" + "  " + "X-Coordinate" + "  " + "Y-Coordinate" + "  " + "Z-Coordinate" + "  "
                        + no.ToString() + "\r\n" + str;
                }
                else
                {
                    str = null + str;
                }
            }
            return str;
        }
        //-------------------------------------------

    }
}
