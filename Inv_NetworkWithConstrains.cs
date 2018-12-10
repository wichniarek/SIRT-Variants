using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc: Network Inversion with Constrains
********************************************/

namespace ClassLibrary_TomoGo
{
    public class Inv_NetworkWithConstrains
    {
        public Inv_NetworkWithConstrains(float[] domain, Inv_SRInfoTable table, List<double> d, InversionParameters inputIP)
        {
            _domain = new float[4];
            _domain[0] = domain[0]; _domain[1] = domain[1]; _domain[2] = domain[2]; _domain[3] = domain[3];
            _nx = inputIP.NX;
            _nz = inputIP.NZ;
            _table = table.copy();

            _iP = new InversionParameters(inputIP.NX, inputIP.NZ);
            _iP.CopyFrom(inputIP);
            conMinX = Vector<double>.Build.Dense(_iP.NX * _iP.NZ, double.MinValue);
            conMaxX = Vector<double>.Build.Dense(_iP.NX * _iP.NZ, double.MaxValue);
            CreateConValues(_iP);// build up conMinValue and conMaxValue
            
            _diffusivity = d.ToList();
            _tMeasured = TsGenerator1(_table);
            _tProcessed = TsGenerator2(_table);//已经乘以系数并开方，也就是说是b项//the vector b
            //MessageBox.Show(_table.F_alpha_d.ToString());

            _fRs = DsGenerator(_domain, _nx, _nz);
            _nodesInv = NodesGenerator8();
            _sirt_Res = new List<SIRT_Result>();
        }
        // Parameter definition
        #region
        //-----------------------------------------------------------
        private InversionParameters _iP;
        private Vector<double> conMinX;
        private Vector<double> conMaxX;
        private float[] _domain;
        private int _nx;
        private int _nz;
        private List<double> _tMeasured;
        private List<double> _tProcessed;
        //-----------------------------------------------------------
        private Inv_SRInfoTable _table;
        private List<double> _diffusivity;
        private List<FRechtangle> _fRs;//通过函数DsGenerator中被赋值//assigned by DsGenerator
        //-----------------------------------------------------------
        private List<NodesInversion> _nodesInv;//通过函数DsGenerator中被赋值//assigned by DsGenerator
        public List<NodesInversion> NodesInv
        {
            get { return _nodesInv; }
            set { _nodesInv = value; }
        }
        //-----------------------------------------------------------
        private List<SIRT_Result> _sirt_Res;
        public List<SIRT_Result> Sirt_Res
        {
            get { return _sirt_Res; }
            set { _sirt_Res = value; }
        }

        public StringBuilder _sbDij_0 = new StringBuilder();
        public StringBuilder _sbDij_0b = new StringBuilder();
        public StringBuilder _sbDij_1b = new StringBuilder();
        public StringBuilder _sbNodes = new StringBuilder();
        public StringBuilder _sbSIRT_1 = new StringBuilder();
        public StringBuilder _sbSIRT_2 = new StringBuilder();
        #endregion
        //--------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------
        ///Measured TravelTime
        ///配置TravelTime的向量-->List<double> TTL//set the vector of travel time
        public List<double> TsGenerator1(Inv_SRInfoTable table)
        {
            List<double> TTL = new List<double>();//--> 经过计算后将用来矩阵计算的TravelTimeList// the vector b used for b=Ax
            foreach (Inv_SRInfoLine L in table.SRInfoLineList)
            {
                Inv_SRInfoLine line = L.copy();
                line.TravelTime = line.TravelTime;
                TTL.Add(line.TravelTime);
            }
            return TTL;
        }
        public List<double> TsGenerator2(Inv_SRInfoTable table)
        {
            List<double> TTL = new List<double>();//--> 经过计算后将用来矩阵计算的TravelTimeList// the vector b used for b=Ax
            foreach (Inv_SRInfoLine L in table.SRInfoLineList)
            {
                Inv_SRInfoLine line = L.copy();
                line.TravelTime = (float)Math.Pow(table.DimFactor * _table.F_alpha_d * line.TravelTime, 0.5);
                TTL.Add(line.TravelTime);
            }
            return TTL;
        }
        //================================================================
        //the following function generates diffusivity rechtangles(Elements)
        public List<FRechtangle> DsGenerator(float[] domain, int nx, int nz)
        {
            List<FRechtangle> List = new List<FRechtangle>();
            if (domain.GetLength(0) != 4 || nx < 1 || nz < 1)
            { MessageBox.Show("Parameter error in DsGenerator"); }
            else
            {
                int elementcount = nx * nz;
                float elementlength = (domain[1] - domain[0]) / nx;
                float elementhigh = (domain[3] - domain[2]) / nz;
                //MessageBox.Show(elementlength.ToString() + "\r\n" + elementhigh.ToString());
                //string str=null;
                for (int i = 0; i < elementcount; i++)
                {
                    int row = i / nx;//show the row where D(i) in;
                    int col = i % nx;//show the column where D(i) in;

                    FRechtangle fr = new FRechtangle(domain[0] + col * elementlength, domain[0] + (col + 1) * elementlength,
                                                     domain[3] - (row + 1) * elementhigh, domain[3] - row * elementhigh);
                    fr.Name = i.ToString();
                    List.Add(fr);
                }
                //MessageBox.Show(str);

            }
            return List;
        }
        //================================================================
        //下面的程序是用来将Grid后的数据类型转换成Dijkstra算法需要的数据类型;
        //the following function is utilized to transfer Grid data into the data format, which is suitable for Dijkstra
        //8 nodes in each cell
        public List<NodesInversion> NodesGenerator8()
        {
            //nx is the number of cells in x direction
            //nz is the number of cells in z direction
            
            _sbNodes.Clear();

            List<NodesInversion> List = new List<NodesInversion>();
            //MessageBox.Show(domain[0].ToString()+"  "+domain[1].ToString()+"  "+domain[2].ToString()+"  "+domain[3].ToString()+"  ");
            if (_domain.GetLength(0) != 4 || _nx < 1 || _nz < 1)
            { MessageBox.Show("Parameter error in NodesGenerator8"); }


            //add common nodes
            int nodescount = 2 * _nx * (_nz + 1) + 2 * _nz * (_nx + 1);
            for (int i = 0; i < nodescount; i++)
            {
                NodesInversion ninv = new NodesInversion();
                ninv.Name = i.ToString();
                List.Add(ninv);
            }
            //设置横向的nodes//set the nodes in x direction
            #region
            for (int i = 0; i < _fRs.Count; i++)
            {
                //下面是设置左上角的node，它下方相邻的矩形是L中的第i个矩形
                List[i * 2].Type = NodesInversion.NodesType.NormalNode;
                List[i * 2].Coor[0] = _fRs[i].Xleft * (2f / 3f) + _fRs[i].Xright * (1f / 3f);
                List[i * 2].Coor[1] = _fRs[i].Zup;
                List[i * 2].Daround[1, 0] = i;
                List[i * 2].Daround[1, 1] = i;
                //下面是设置右上角的node，它下方相邻的矩形是L中的第i个矩形
                List[i * 2 + 1].Type = NodesInversion.NodesType.NormalNode;
                List[i * 2 + 1].Coor[0] = _fRs[i].Xleft * (1f / 3f) + _fRs[i].Xright * (2f / 3f);
                List[i * 2 + 1].Coor[1] = _fRs[i].Zup;
                List[i * 2 + 1].Daround[1, 0] = i;
                List[i * 2 + 1].Daround[1, 1] = i;
                //下面是设置左下角的node，它上方相邻的矩形是L中的第i个矩形
                List[i * 2 + 2 * _nx].Type = NodesInversion.NodesType.NormalNode;
                List[i * 2 + 2 * _nx].Coor[0] = _fRs[i].Xleft * (2f / 3f) + _fRs[i].Xright * (1f / 3f);
                List[i * 2 + 2 * _nx].Coor[1] = _fRs[i].Zdown;
                List[i * 2 + 2 * _nx].Daround[0, 0] = i;
                List[i * 2 + 2 * _nx].Daround[0, 1] = i;
                //下面是设置右下角的node，它上方相邻的矩形是L中的第i个矩形
                List[i * 2 + 2 * _nx + 1].Type = NodesInversion.NodesType.NormalNode;
                List[i * 2 + 2 * _nx + 1].Coor[0] = _fRs[i].Xleft * (1f / 3f) + _fRs[i].Xright * (2f / 3f);
                List[i * 2 + 2 * _nx + 1].Coor[1] = _fRs[i].Zdown;
                List[i * 2 + 2 * _nx + 1].Daround[0, 0] = i;
                List[i * 2 + 2 * _nx + 1].Daround[0, 1] = i;

            }
            #endregion
            //设置竖向的nodes//set the nodes in z direction
            #region

            int m = 2 * _nx * (_nz + 1);//已用nodes数量，在上面一个region中已经标记的nodes的数量
            int count = 0;//用来计数，表明在下面双重循环的次数
            for (int i = 0; i < _nx; i++)//大循环，下标第二位
            {

                for (int j = 0; j < _nz; j++)//小循环，下标第一位
                {

                    int D = j * _nx + i;//设置D的序号
                    int node0 = m + 2 * count;//左上角第一个node
                    //MessageBox.Show("i="+i.ToString()+" j="+j.ToString()+"  D="+D.ToString()+"  node0="+node0.ToString());
                    //下面是左上角node的设置
                    List[node0].Type = NodesInversion.NodesType.NormalNode;
                    List[node0].Coor[0] = _fRs[D].Xleft;
                    List[node0].Coor[1] = _fRs[D].Zup * (2f / 3f) + _fRs[D].Zdown * (1f / 3f);
                    List[node0].Daround[0, 1] = D;//右边相邻的矩形是L中第D个矩形
                    List[node0].Daround[1, 1] = D;//右边相邻的矩形是L中第D个矩形
                    //下面是右上角node的设置
                    List[node0 + 1].Type = NodesInversion.NodesType.NormalNode;
                    List[node0 + 1].Coor[0] = _fRs[D].Xleft;
                    List[node0 + 1].Coor[1] = _fRs[D].Zup * (1f / 3f) + _fRs[D].Zdown * (2f / 3f);
                    List[node0 + 1].Daround[0, 1] = D;//左边相邻的矩形是L中第D个矩形
                    List[node0 + 1].Daround[1, 1] = D;//左边相邻的矩形是L中第D个矩形
                    //下面是左下角node的设置
                    List[node0 + 2 * _nz].Type = NodesInversion.NodesType.NormalNode;
                    List[node0 + 2 * _nz].Coor[0] = _fRs[D].Xright;
                    List[node0 + 2 * _nz].Coor[1] = _fRs[D].Zup * (2f / 3f) + _fRs[D].Zdown * (1f / 3f);
                    List[node0 + 2 * _nz].Daround[0, 0] = D;//右边相邻的矩形是L中第D个矩形
                    List[node0 + 2 * _nz].Daround[1, 0] = D;//右边相邻的矩形是L中第D个矩形
                    //下面是右下角node的设置
                    List[node0 + 2 * _nz + 1].Type = NodesInversion.NodesType.NormalNode;
                    List[node0 + 2 * _nz + 1].Coor[0] = _fRs[D].Xright;
                    List[node0 + 2 * _nz + 1].Coor[1] = _fRs[D].Zup * (1f / 3f) + _fRs[D].Zdown * (2f / 3f);
                    List[node0 + 2 * _nz + 1].Daround[0, 0] = D;//左边相邻的矩形是L中第D个矩形
                    List[node0 + 2 * _nz + 1].Daround[1, 0] = D;//左边相邻的矩形是L中第D个矩形

                    count = count + 1;
                }
            }

            #endregion
            // set S
            #region
            List<NodesInversion> List_temp = new List<NodesInversion>();
            //MessageBox.Show(table.Sset.Count.ToString());
            foreach (Inv_SR3D sender in _table.Sset)
            {
                bool criterium = false;
                foreach (NodesInversion sender2 in List)
                {

                    if (sender.Coor[0] == sender2.Coor[0] && sender.Coor[2] == sender2.Coor[1])
                    {
                        sender2.Name = sender2.Name + "(" + sender.Name + ")"; break;
                    }
                    else
                    {
                        criterium = true;
                    }
                }
                if (criterium)
                {
                    NodesInversion newnode = new NodesInversion();
                    newnode.Name = sender.Name;
                    newnode.Type = NodesInversion.NodesType.Sender;
                    newnode.Coor[0] = sender.Coor[0];
                    newnode.Coor[1] = sender.Coor[2];
                    for (int i = 0; i < _fRs.Count; i++)
                    {

                        int k = _fRs[i].Location(newnode.Coor[0], newnode.Coor[1]);

                        switch (k)
                        {
                            case -1: MessageBox.Show("error by using rechtsangle.loction"); break;
                            case 0: newnode.Daround[1, 1] = i; break;
                            case 1: newnode.Daround[1, 0] = i; newnode.Daround[1, 1] = i; break;
                            case 2: newnode.Daround[1, 0] = i; break;
                            case 3: newnode.Daround[0, 1] = i; newnode.Daround[1, 1] = i; break;
                            case 4: newnode.Daround[0, 0] = i; newnode.Daround[0, 1] = i; newnode.Daround[1, 0] = i; newnode.Daround[1, 1] = i; break;
                            case 5: newnode.Daround[0, 0] = i; newnode.Daround[1, 0] = i; break;
                            case 6: newnode.Daround[0, 1] = i; break;
                            case 7: newnode.Daround[0, 0] = i; newnode.Daround[0, 1] = i; break;
                            case 8: newnode.Daround[0, 0] = i; break;
                            case 9: break;
                        }

                    }
                    List_temp.Add(newnode);
                }

            }

            #endregion
            //set R
            #region
            foreach (Inv_SR3D rec in _table.Rset)
            {
                bool criterium = false;
                foreach (NodesInversion rec2 in List)
                {
                    if (rec.Coor[0] == rec2.Coor[0] && rec.Coor[2] == rec2.Coor[1])
                    {
                        rec2.Name = rec2.Name + "(" + rec.Name + ")"; break;
                    }
                    else
                    {
                        criterium = true;
                    }
                }
                if (criterium)
                {
                    NodesInversion newnode = new NodesInversion();
                    newnode.Name = rec.Name;
                    newnode.Type = NodesInversion.NodesType.Receiver;
                    newnode.Coor[0] = rec.Coor[0];
                    newnode.Coor[1] = rec.Coor[2];
                    for (int i = 0; i < _fRs.Count; i++)
                    {
                        //MessageBox.Show(i.ToString());
                        int k = _fRs[i].Location(newnode.Coor[0], newnode.Coor[1]);
                        switch (k)
                        {
                            case -1: MessageBox.Show("error by using rechtsangle.loction"); break;
                            case 0: newnode.Daround[1, 1] = i; break;
                            case 1: newnode.Daround[1, 0] = i; newnode.Daround[1, 1] = i; break;
                            case 2: newnode.Daround[1, 0] = i; break;
                            case 3: newnode.Daround[0, 1] = i; newnode.Daround[1, 1] = i; break;
                            case 4: newnode.Daround[0, 0] = i; newnode.Daround[0, 1] = i; newnode.Daround[1, 0] = i; newnode.Daround[1, 1] = i; break;
                            case 5: newnode.Daround[0, 0] = i; newnode.Daround[1, 0] = i; break;
                            case 6: newnode.Daround[0, 1] = i; break;
                            case 7: newnode.Daround[0, 0] = i; newnode.Daround[0, 1] = i; break;
                            case 8: newnode.Daround[0, 0] = i; break;
                            case 9: break;
                        }
                    }
                    List_temp.Add(newnode);
                }
            }
            #endregion
            //copy S and R into the List
            foreach (NodesInversion n in List_temp)
                List.Add(n);
            //
            foreach (NodesInversion t in List)
            {
                _sbNodes.Append(t.print(1));
            }
            return List;
        }
        //================================================================
        //the following function is for Dijkstra algorithm
        //add nodes
        public List<Edge> Dij_0(List<NodesInversion> input)
        {
            _sbDij_0.Clear();

            //NodesInversion to Node
            Dictionary<string, Node> dictnode = new Dictionary<string, Node>();
            foreach (NodesInversion n in input)
            {
                Node newnode = new Node(n.Name);
                dictnode.Add(n.Name, newnode);
            }


            //define the output
            List<Edge> output = new List<Edge>();
            //设置“有距离”点之间的距离
            //set the distance between nodes whose distance is large than 0
            foreach (NodesInversion n in input)
            {
                foreach (NodesInversion m in input)
                {
                    bool criterium = false;
                    if (!string.Equals(n.Name, m.Name))//n!=m. n and m are different nodes
                    {
                        List<int> nDraund = new List<int>();
                        nDraund.Add(n.Daround[0, 0]); nDraund.Add(n.Daround[0, 1]); nDraund.Add(n.Daround[1, 0]); nDraund.Add(n.Daround[1, 1]);
                        List<int> mDraund = new List<int>();
                        mDraund.Add(m.Daround[0, 0]); mDraund.Add(m.Daround[0, 1]); mDraund.Add(m.Daround[1, 0]); mDraund.Add(m.Daround[1, 1]);
                        foreach (int i in nDraund)
                        {
                            if (mDraund.Contains(i) && i != -1)
                            {
                                criterium = true;
                                break;
                            }
                        }
                        if (nDraund.SequenceEqual(mDraund))
                            criterium = false;
                    }
                    if (criterium)
                    {
                        double distance = Math.Pow(Math.Pow(n.Coor[0] - m.Coor[0], 2) + Math.Pow(n.Coor[1] - m.Coor[1], 2), 0.5);
                        Edge newedge = new Edge(dictnode[n.Name], dictnode[m.Name], distance);
                        output.Add(newedge);
                    }
                }

            }

            // record process
            foreach (Edge t in output)
            {
                _sbDij_0.Append(t.print(1));
            }
            
            return output;
        }
        //the following function is based on the previous function, the D will be added into Edge
        public List<Edge> Dij_0b(List<NodesInversion> input_Nodes, List<double> input_D)
        {

            _sbDij_0b.Clear();
            _sbDij_0b.Append("*** log file of function Dij_0b ***");
            _sbDij_0b.Append("----------------------------\r\n");
            //================================================================
            //Calculate the root of D
            Dictionary<int, double> velocity = new Dictionary<int, double>();
            if (input_D.Count != 0)
                for (int i = 0; i < input_D.Count; i++)
                {
                    velocity.Add(i, 1.0 / (Math.Pow(input_D[i], 0.5)));
                }
            else
                MessageBox.Show("parameter error in Dij_0b");

            
            _sbDij_0b.Append("input   velocity\r\n");
            for (int i = 0; i < input_D.Count; i++)
            {
                _sbDij_0b.Append(input_D[i].ToString() + " " + velocity[i].ToString() + "\r\n");
            }
            
            //================================================================d
            //NodesInversion to Node
            Dictionary<string, Node> dictnode = new Dictionary<string, Node>();
            foreach (NodesInversion n in input_Nodes)
            {
                Node newnode = new Node(n.Name);
                dictnode.Add(n.Name, newnode);
            }


            //define output
            List<Edge> output = new List<Edge>();
            //设置“有距离”点之间的距离
            //set the distance between nodes whose distance is large than 0
            foreach (NodesInversion n in input_Nodes)
            {
                foreach (NodesInversion m in input_Nodes)
                {
                    bool criterion = false;
                    int markD = -1;
                    if (!string.Equals(n.Name, m.Name))//n!=m. n and m are different nodes
                    {
                        List<int> nA = new List<int>();
                        nA.Add(n.Daround[0, 0]); nA.Add(n.Daround[0, 1]); nA.Add(n.Daround[1, 0]); nA.Add(n.Daround[1, 1]);
                        List<int> mA = new List<int>();
                        mA.Add(m.Daround[0, 0]); mA.Add(m.Daround[0, 1]); mA.Add(m.Daround[1, 0]); mA.Add(m.Daround[1, 1]);

                        if (!nA.SequenceEqual(mA))
                        {
                            foreach (int i in nA)
                            {
                                if (mA.Contains(i) && i != -1)
                                {
                                    criterion = true;
                                    markD = i;
                                    break;
                                }
                            }

                        }
                    }
                    if (criterion)
                    {
                        double distance = Math.Pow(Math.Pow((double)n.Coor[0] - (double)m.Coor[0], 2.0) + Math.Pow((double)n.Coor[1] - (double)m.Coor[1], 2.0), 0.5);
                        //MessageBox.Show(distance.ToString());
                        Edge newedge = new Edge(dictnode[n.Name], dictnode[m.Name], distance * velocity[markD]);
                        output.Add(newedge);
                    }
                }

            }

            //print to txt
            if (output.Count > 0)
                foreach (Edge t in output)
                {
                    _sbDij_0b.Append(t.print(1));
                }
            //System.IO.File.WriteAllText("Inv_NetworkWithConstrains--Dij_0b.txt", _sbDij_0b.ToString());
            return output;
        }
        
        //下面的程序适用Dijkstra中的S到R情况,Dij_1b调用了Dij_0b
        //the following function is for the Dijkstra S to R, Dij_1b uses Dij_0b
        public List<Dijkstra_output> Dij_1b(List<NodesInversion> input, List<double> D, Inv_SRInfoTable inputTable)
        {
            _sbDij_1b.Clear();
            //==============================================================
            //定义输出//define output
            List<Dijkstra_output> output = new List<Dijkstra_output>();
            //将乘以D后变成以根号6ftor4ft为单位的距离输出到Edge中
            //D*root(6ftor4ft) as distance and write into Edge
            List<Edge> _edges_t = Dij_0b(input, D);

            //==============================================================
            Dictionary<string, Node> _dictNodes = new Dictionary<string, Node>();
            List<Node> _nodes = new List<Node>();

            foreach (NodesInversion n in input)
            {
                _nodes.Add(new Node(n.Name));
                _dictNodes.Add(n.Name, new Node(n.Name));
            }
            //将S和R分别放在ListS和ListR中
            //write S and R into ListS and ListR, perspectively
            List<Node> ListS = new List<Node>();
            List<Node> ListR = new List<Node>();
            #region
            foreach (NodesInversion n in input)
            {
                switch (n.Type)
                {
                    case NodesInversion.NodesType.Sender:
                        ListS.Add(new Node(n.Name)); break;
                    case NodesInversion.NodesType.Receiver:
                        ListR.Add(new Node(n.Name)); break;
                    default:
                        break;
                }
            }
            #endregion

            foreach (Node s in ListS)
            {
                // new Dijkstra class 
                Dijkstra d = new Dijkstra(_edges_t, _nodes);
                // calculate the distance
                d.calculateDistance(s);
                foreach (Node r in ListR)
                {                    
                    // here need to check if s->r we need???
                    bool sr = false;
                    foreach(Inv_SRInfoLine line in inputTable.SRInfoLineList)
                    {
                        if(line.Start.Name.Equals(s.Name)&&line.End.Name.Equals(r.Name))
                        {
                            sr = true;
                            break;
                        }
                    }
                    if (sr)
                    {
	                    // Pfad zu einem bestimmten Knoten ausgeben
	                    List<Node> path = d.getPathTo(r);
	                    Dijkstra_output dij = new Dijkstra_output();
	                    dij.Start = s;
	                    dij.End = r;
	                    dij.Path = path;
	                    dij.TravelTime = d.Dist[r.Name];
	                    output.Add(dij);
                    }
                }
            }
            // log file 
            foreach (Dijkstra_output n in output)
            {
                _sbDij_1b.Append(n.print());
            }

            return output;

        }
        //下面是SIRT主程序//the main function of SIRT
        public void SIRT_2()
        {
            _sbSIRT_2.Clear();
            //========================================================================
            List<Edge> sirt_edges = Dij_0(_nodesInv);//将以距离作为单位的edge包装好，准备调用//prepair the edge
            //========================================================================
            Dictionary<string, NodesInversion> sirt_nodesinv = new Dictionary<string, NodesInversion>();
            foreach (NodesInversion n in _nodesInv)//将nodesInv包装好，适合通过name调用//prepair the nodesInv, for usage by name
                sirt_nodesinv.Add(n.Name, n);
            //========================================================================
            int n_d = _diffusivity.Count;//距离矩阵A的列数，D的数量//the number of columns in A, equals the number of D
            int n_t = _tProcessed.Count;//距离矩阵A的行数，T的数量//the number of rows in A, equals the number of T
            #region
            _sbSIRT_2.Append("TomoGo 1.01" + "\r\n");
            _sbSIRT_2.Append("f: " + _table.F_alpha_d.ToString() + "\r\n");
            _sbSIRT_2.Append("count of D (col of distance-matrix A): " + n_d.ToString() + "\r\n");
            _sbSIRT_2.Append("count of T (row of distance-matrix A): " + n_t.ToString() + "\r\n");
            _sbSIRT_2.Append("Initial Diffusivity:\r\n");
            foreach (double d in _diffusivity)
            {
                _sbSIRT_2.Append(d.ToString("N") + "    ");
            }
            _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
            _sbSIRT_2.Append("Initial x=1/sqrt(D):\r\n");
            foreach (double d in _diffusivity)
            {
                _sbSIRT_2.Append((1/Math.Sqrt(d)).ToString("N") + "    ");
            }
            _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
            _sbSIRT_2.Append("Measured Travel Time t:\r\n");
            foreach (double t in _tMeasured)
            {
                _sbSIRT_2.Append(t.ToString("N") + "    ");
            }
            _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
            _sbSIRT_2.Append("Measured b = sqrt(6*f*t or 4*f*t):\r\n");
            foreach (double t in _tProcessed)
            {
                _sbSIRT_2.Append(t.ToString("N") + "    ");
            }
            _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
            #endregion
            //====================================================================
            // initialize b = sqrt(6*f*t or4*f*t)
            Vector<double> b = Vector<double>.Build.DenseOfArray(_tProcessed.ToArray());
            // initialize x, x=1/(root of diffusivity)
            Vector<double> x = Vector<double>.Build.Dense(_diffusivity.Count, 1);
            for (int i = 0; i < _diffusivity.Count; i++)
            {
                x[i] = 1.0 / (Math.Pow(_diffusivity[i], 0.5));
            }
            List<double> diffusivity = _diffusivity.ToList();
            //================================================================
            #region
            int iter = 0;
            double RMS = _iP.CriRMS + 1;
            while (iter < _iP.CurRay & RMS > _iP.CriRMS)
            {
                SIRT_Result si_R = new SIRT_Result();
                si_R.Iter = iter;//将用时写入output中//wirte traveltime into output
                //================================================================
                si_R.GuessDiffusivity = diffusivity.ToList();//将Guess的D写入output中//write D of Guess into output
                //================================================================
                #region
                _sbSIRT_2.Append("************************************" + "\r\n");
                _sbSIRT_2.Append("**********  iteration: " + iter.ToString() + "  **********" + "\r\n");
                _sbSIRT_2.Append("************************************" + "\r\n");
                _sbSIRT_2.Append("Input Diffusivity: " + "\r\n");
                foreach (double d in diffusivity)
                {
                    _sbSIRT_2.Append(d.ToString("N") + "    ");
                }
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                //===============================================================
                _sbSIRT_2.Append("************************************" + "\r\n");
                _sbSIRT_2.Append("**********  iteration: " + iter.ToString() + "  **********" + "\r\n");
                _sbSIRT_2.Append("************************************" + "\r\n");
                _sbSIRT_2.Append("d_cori: " + "\r\n");
                foreach (double d in diffusivity)
                {
                    _sbSIRT_2.Append(d.ToString("N") + "    ");
                }
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                #endregion
                //===============================================================
                //STEP 0
                List<Dijkstra_output> dij_outputs = Dij_1b(_nodesInv, diffusivity,_table);//每次迭代产生的Dijkstra算法的结果//the result of Dijkstra algorithm after every iteration
                si_R.Dijk_output = dij_outputs.ToList();//将用时写入output中//write traveltime into output
                #region
                foreach (Dijkstra_output n in dij_outputs)
                {
                    _sbSIRT_2.Append(n.print());
                }
                #endregion
                //====================================================
                //STEP 1
                //通过对Dijkstra Algorithm结果的逐行分析，设置距离矩阵。
                //analyze the Dijkstra algorithm result, set the distance matrix A
                Matrix<double> A = Matrix<double>.Build.Dense(n_t, n_d, 0);
                #region
                int row = 0;
                foreach (Dijkstra_output dij_single in dij_outputs)//对Dij算法的结果逐行分析//Dijkstra result analysis
                {
                    List<Node> path = dij_single.Path;

                    if (path.Count < 2)
                        MessageBox.Show("error in path.count in SIRT_2");
                    for (int i = 0; i < path.Count - 1; i++)
                    {
                        Node n0 = path[i];
                        Node n1 = path[i + 1];
                        NodesInversion n0_inv = sirt_nodesinv[n0.Name];
                        NodesInversion n1_inv = sirt_nodesinv[n1.Name];
                        List<int> union = n0_inv.union(n1_inv);
                        if (union != null)
                        {
                            foreach (Edge e in sirt_edges)
                            {
                                if (e.Origin.equals(n0) && e.Destination.equals(n1))
                                {
                                    A[row, union[0]] = e.Distance;
                                    break;
                                }
                            }
                        }
                    }
                    row = row + 1;
                }
                #endregion
                _sbSIRT_2.Append("***** Ax=b, x->x2*****\r\n");
                _sbSIRT_2.Append("A: \r\n" + A.ToString() + "\r\n------------------------------------\r\n");
                //====================================================
                _sbSIRT_2.Append("x: \r\n" + x.ToString() + "\r\n------------------------------------\r\n");
                //====================================================
                //STEP 3 : calculate c=Ax
                //计算用时//calculate the traveltime
                Vector<double> c = A * x;
                #region
                _sbSIRT_2.Append("c=A*x: " + "\r\n");
                _sbSIRT_2.Append(c.ToString("N") + "   ");
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                #endregion
                //====================================================
                //TEST:
                //从Dijkstra中提取时间, compare the time with Step 4 where the time is from A*x
                //get traveltime from Dijkstra, compare the time with Step 4 where the time is from A*x
                double[] t_dijkstra = new double[n_t];
                for (int i = 0; i < n_t; i++)
                    t_dijkstra[i] = dij_outputs[i].TravelTime;
                #region
                _sbSIRT_2.Append("t_dijkstra: " + "\r\n");
                foreach (double dou in t_dijkstra)
                {
                    _sbSIRT_2.Append(dou.ToString("N") + "   ");
                }
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                _sbSIRT_2.Append("real time from dijkstra: " + "\r\n");
                foreach (double dou in t_dijkstra)
                {
                    _sbSIRT_2.Append((dou * dou / (_table.F_alpha_d * _table.DimFactor)).ToString("N") + "   ");
                }
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                #endregion

                //====================================================
                //STEP 4
                //计算时间差//calculate the delta t
                double[] realTimeList = new double[n_t]; //转换时间为真实的时间  //transfer to real time
                for (int i = 0; i < n_t; i++)
                {
                    double realTime = Math.Pow(c[i], 2) / ((double)_table.DimFactor * (double)_table.F_alpha_d);//the realtime
                    realTimeList[i] = realTime;
                    si_R.Time.Add(realTime);//将用时写入output中//wirte traveltime into output
                }
                Statistic stc = new Statistic();
                RMS = stc.RMS(new List<double>(realTimeList), _tMeasured);
                #region
                _sbSIRT_2.Append("Real Time from c=A*x: " + "\r\n");
                foreach (double rt in realTimeList)
                {
                    _sbSIRT_2.Append(rt.ToString("N") + "   ");
                }
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                _sbSIRT_2.Append("RMS Residual: " + RMS.ToString("N"));
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                #endregion
                //====================================================
                //STEP 6
                SIRT_Options so = new SIRT_Options(b, A, x, 1);
                Vector<double> x2;
                switch(_iP.SIRTOption)
                {
                    case 1:
                        x2 = so.XUpdaterDROP(x, 1);//use the DROP method
                        break;
                    case 2:
                        x2 = so.XUpdaterCAV(x, 1);//use the CAV method
                        break;
                    case 3:
                        x2 = so.XUpdaterCimmino(x,1); //use the Cimmino method              
                        break; 
                    default:
                        x2 = x.Clone();// no update
                        break;
                }
                #region
                _sbSIRT_2.Append("x2: " + "\r\n");
                _sbSIRT_2.Append(x2.ToString("N") + "   ");
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                #endregion
                //====================================================
                //STEP 7: update x and diffusivity
                x = x2.Clone();
                for (int i = 0; i < n_d; i++)
                {
                    if (x[i] < conMinX[i])
                        x[i] = conMinX[i];
                    if (x[i] > conMaxX[i])
                        x[i] = conMaxX[i];
                    diffusivity[i] = 1.0 / (x[i] * x[i]);
                }
                #region
                _sbSIRT_2.Append("x after Min&Max filter: " + "\r\n");
                _sbSIRT_2.Append(x.ToString("N") + "   ");
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                _sbSIRT_2.Append("diffusivity after Min&Max filter: " + "\r\n");
                foreach(double dou in diffusivity)
                    _sbSIRT_2.Append(dou.ToString("N") + "   ");
                _sbSIRT_2.Append("\r\n" + "------------------------------------" + "\r\n");
                #endregion
                //====================================================
                iter = iter + 1;
                //====================================================
                _sbSIRT_2.Append("====================================\r\n");
                si_R.UpdatedDiffusivity = diffusivity.ToList();//write the result into output
                _sirt_Res.Add(si_R);//write the result into output

            }
            #endregion
            //================================================================
            // iteration is finished
            _diffusivity = diffusivity.ToList();//update the diffusivity
            //to check the result in txt file, see below
            //System.IO.File.WriteAllText("Inversion-Network-SIRT2.txt", sbPrint.ToString());
            //System.IO.File.WriteAllText("Inversion-Network-SIRT2-Dijk.txt", _sbSIRT_2.ToString());

        }
        private void CreateConValues(InversionParameters inputIP)
        {
            double xMax = 1.0 / (Math.Sqrt(_iP.MinD));// set the constrain of x, x=1/(root of diffusivity)
            double xMin = 1.0 / (Math.Sqrt(_iP.MaxD));
            double xBig = xMax;
            double xSma = xMin;
            int index = 0;
            int rows = inputIP.ConBoolMatrix.RowCount;
            int cols = inputIP.ConBoolMatrix.ColumnCount;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    index = i * cols + j;
                    // every cell has constrin: >xMin and <xMax
                    conMinX[index] = xMin;
                    conMaxX[index] = xMax;
                    if (inputIP.ConBoolMatrix[i, j] == 1)
                    {
                        // this cell has a extra constrain
                        xBig = 1.0 / (Math.Sqrt(_iP.ConValueMatrix[i, j]));
                        xSma = 1.0 / (Math.Sqrt(_iP.ConValueMatrix[i, j]));
                        conMinX[index] = xSma;
                        conMaxX[index] = xBig;
                    }
                    else if (inputIP.ConBoolMatrix[i, j] == 0)
                    {
                        // this cell does not any constrain more  
                    }
                    else
                    { MessageBox.Show("Error in CreateConValues in Inv_StraightRayWithConstrins"); }
                }
            }
            // check if every value in conMaxX is larger than conMinX
            for (int i = 0; i < conMaxX.Count; i++)
            {
                if (conMaxX[i] < conMinX[i])
                {
                    MessageBox.Show("Error in Check of CreateConValues in Inv_StraightRayWithConstrins");
                }
            }

        }
    }
}
