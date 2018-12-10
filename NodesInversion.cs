using System;
using System.Collections.Generic;
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
    public class NodesInversion
    {
        //********************************
        public NodesInversion()
        {
            _name = null;// the name of the node
            _type = NodesType.Undefined; // differnent kinds of node
            _coor[0] = float.MaxValue; _coor[1] = float.MaxValue;// the coordinate of node
            _daround[0, 0] = -1; _daround[0, 1] = -1; _daround[1, 0] = -1; _daround[1, 1] = -1; // the four quadrants of the node
        }
        //-------------------------------
        private string _name;
        public string Name
        {
            get { return _name; }
            set { _name = value; }
        }
        //-------------------------------
        public enum NodesType
        {
            Receiver, Sender, NormalNode, Undefined
        }
        private NodesType _type;
        public NodesType Type
        {
            get { return _type; }
            set { _type = value; }
        }
        //-------------------------------
        private float[] _coor = new float[2];
        public float[] Coor
        {
            get { return _coor; }
            set { _coor = Coor; }
        }
        //-------------------------------
        private int[,] _daround = new int[2, 2];
        public int[,] Daround
        {
            get { return _daround; }
            set { _daround = value; }
        }
        //-------------------------------

        //-------------------------------
        //下面的函数用来计算m和n之间的连线所在的Element
        //the following function is used to identify the cell where the node m, n and their connection stay.
        public List<int> union(NodesInversion n)
        {
            List<int> output = new List<int>();
            //限定：连线只能处于一个Element中
            //constraint: the connection between nodes can only be in one cell
            List<int> Draund = new List<int>();
            Draund.Add(_daround[0, 0]); Draund.Add(_daround[0, 1]); Draund.Add(_daround[1, 0]); Draund.Add(_daround[1, 1]);
            List<int> nDraund = new List<int>();
            nDraund.Add(n.Daround[0, 0]); nDraund.Add(n.Daround[0, 1]); nDraund.Add(n.Daround[1, 0]); nDraund.Add(n.Daround[1, 1]);
            foreach (int i in Draund)
            {
                foreach (int j in nDraund)
                {
                    if (i == j && i != -1)
                    {
                        output.Add(i);
                        nDraund.Remove(j);
                        break;
                    }

                }
            }
            return output;//output=null表示没有union的element// output=null means there isn't any joint cell
        }
        //-------------------------------

        //-------------------------------
        public string print(int n)
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
                str = _name + " "
                    + _type.ToString() + "    "
                    + _coor[0].ToString() + " " + _coor[1].ToString() + " "
                    + _daround[0, 0].ToString() + " " + _daround[0, 1].ToString() + " " + _daround[1, 0].ToString() + " " + _daround[1, 1].ToString() + "\r\n";
                if (n == 0)
                {
                    str = "Name" + "    " + "Type" + "  " + "X-Coordinate" + "  " + "Y-Coordinate" + "  "
                        + "D_00" + "   " + "D_01" + "   " + "D_10" + "   " + "D_11" + "   " + "\r\n" + str;
                }
                else
                {
                    str = null + str;
                }
            }
            return str;
        }
    }
}
