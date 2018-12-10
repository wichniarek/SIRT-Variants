using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc: The Input Table
********************************************/

namespace ClassLibrary_TomoGo
{
    public class Inv_SRInfoTable
    {
        //****************************************
        public Inv_SRInfoTable()
        {
            sRInfoLineList = new ArrayList();
            _dimFactor = 6; //6 for 3D,4 for 2D 
            _p = 100;//100 for t100
            f_alpha_d = 1; //for TPeak
            sset = new List<Inv_SR3D>();
            rset = new List<Inv_SR3D>();
        }
        //****************************************
        private ArrayList sRInfoLineList;
        public ArrayList SRInfoLineList
        {
            get { return sRInfoLineList; }
            set { sRInfoLineList = value; }
        }
        //-----------------------------------------  
        private float _dimFactor;
        public float DimFactor
        {
            get { return _dimFactor; }
            set { _dimFactor = value; }
        }
        //-----------------------------------------
        private int _p;
        public int P
        {
            get { return _p; }
            set { _p = value; }
        }
        //-----------------------------------------
        private float f_alpha_d;
        public float F_alpha_d
        {
            get { return f_alpha_d; }
            set { f_alpha_d = value; }
        }
        //-----------------------------------------    
        private List<Inv_SR3D> sset;
        public List<Inv_SR3D> Sset
        {
            get { return sset; }
            set { sset = value; }
        }
        //-----------------------------------------   
        private List<Inv_SR3D> rset;
        public List<Inv_SR3D> Rset
        {
            get { return rset; }
            set { rset = value; }
        }
        //-----------------------------------------    

        //----------------------------------------- 
        //****************************************
        public Inv_SRInfoTable copy()
        {
            Inv_SRInfoTable output = new Inv_SRInfoTable();
            foreach (Inv_SRInfoLine line in this.sRInfoLineList)
                output.sRInfoLineList.Add(line);
            output._dimFactor = this._dimFactor;
            output.f_alpha_d = this.f_alpha_d;
            output._p = this._p;
            output.sset = this.sset.ToList();
            output.rset = this.rset.ToList();
            return output;
        }
        //----------------------------------------- 
        public int countt()
        {
            //countt is a function, 
            //the return value is the number of traveltimes recorded in table
            int n = SRInfoLineList.Count;
            return n;
        }
        public int counts()
        {
            //counts is a function,
            //the return value is the number of senders recorded in table
            int n = 0;
            ArrayList list = new ArrayList();
            foreach (Inv_SRInfoLine line in SRInfoLineList)
            {
                String name = line.Start.Name;
                if (list.Contains(name) == false)
                {
                    n = n + 1;
                    list.Add(name);
                }
            }
            return n;
        }
        public int countr()
        {
            //counts is a function,
            //the return value is the number of receivers recorded in table
            int n = 0;
            ArrayList list = new ArrayList();
            foreach (Inv_SRInfoLine line in SRInfoLineList)
            {
                String name = line.End.Name;
                if (list.Contains(name) == false)
                {
                    n = n + 1;
                    list.Add(name);
                }
            }
            return n;
        }
        public float[] domain()
        {
            float[] array = new float[4];

            float x1, x2, z1, z2;//x1<x2;z1<z2;
            x1 = ((Inv_SRInfoLine)(SRInfoLineList[0])).Start.Coor[0]; x2 = x1;
            z1 = ((Inv_SRInfoLine)(SRInfoLineList[0])).Start.Coor[2]; z2 = z1;
            int n = SRInfoLineList.Count;
            for (int i = 0; i < n; i++)
            {
                float x1t = ((Inv_SRInfoLine)(SRInfoLineList[i])).Start.Coor[0];
                float x2t = ((Inv_SRInfoLine)(SRInfoLineList[i])).End.Coor[0];
                float z1t = ((Inv_SRInfoLine)(SRInfoLineList[i])).Start.Coor[2];
                float z2t = ((Inv_SRInfoLine)(SRInfoLineList[i])).End.Coor[2];
                if (x1 > x1t) x1 = x1t;
                if (x1 > x2t) x1 = x2t;
                if (z1 > z1t) z1 = z1t;
                if (z1 > z2t) z1 = z2t;
                if (x2 < x1t) x2 = x1t;
                if (x2 < x2t) x2 = x2t;
                if (z2 < z1t) z2 = z1t;
                if (z2 < z2t) z2 = z2t;
            }

            array[0] = x1;
            array[1] = x2;
            array[2] = z1;
            array[3] = z2;

            return array;
        }
        //----------------------------------------- 
        public string print()
        {
            string str = "SRInfoLineList:\r\n";
            foreach (Inv_SRInfoLine m in sRInfoLineList)
            {
                str = str + m.Print(1);
            }
            str = str + "\r\nDimFactor:\r\n";
            str=  str + _dimFactor.ToString() + "  ";
            str = str + "\r\nP:\r\n";
            str = str + _p.ToString() + "  ";
            str = str + "\r\nf_alpha_d:\r\n";
            str = str + f_alpha_d.ToString() + "  ";
            return str;
        }

    }
}
