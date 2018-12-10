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
    public class Main
    {
        public Main()
        {
            _domain[0] = 0;
            _domain[1] = 1;
            _domain[2] = 0;
            _domain[3] = 1;
          
        }
        //--------------------------------------------
        private float[] _domain;//the domain of the investigation area
        public float[] Domain
        {
            get { return _domain; }
            set { _domain = value; }
        }
        private List<double> _d;//the initial diffusivity distribution as vector
        public List<double> D
        {
            get { return _d; }
            set { _d = value; }
        }
        private Inv_SRInfoTable _table;//the table contains information of signals and nodes.
        public Inv_SRInfoTable Table
        {
            get { return _table; }
            set { _table = value; }
        }
        private InversionParameters _iP;//this class contains information of some constraints
        public InversionParameters IP
        {
            get { return _iP; }
            set { _iP = value; }
        }

        //--------------------------------------------
        public void SIRT()
        {
            //===================================================== 
            int nx = _iP.NX;//the resolution in x direction
            int nz = _iP.NZ;//the resolution in z direction
            int A = _iP.StrRay;//the number of straight ray iteration
            int B = _iP.CurRay;//the number of curved ray iteration
            //=====================================================
            //SETP 1 : Straight Ray Algorithm
            List<double> D1 = new List<double>();
            List<double> D2 = new List<double>();
            for (int i = 0; i < _iP.ConValueMatrix.RowCount; i++)
            {
                for (int j = 0; j < _iP.ConValueMatrix.ColumnCount; j++)
                {
                    D1.Add(_iP.ConValueMatrix[i, j]);//copy initial value to D1
                }
            }
            if (A != 0)
            {
                Inv_StraightRayWithConstrain inv_strray = new Inv_StraightRayWithConstrain(_domain, _table, D1, _iP);
                inv_strray.StRayAlgo();// iteration=0 is fine
                _sIRT_Res_StrRay = inv_strray.Result.ToList();
                int iterInvStr = inv_strray.Result.Count;
                D2 = inv_strray.Result[iterInvStr - 1].newD2newD().ToList();//this is result D from Straight Ray, and also is initial D for Newtwork algorithm
                //_sbStrRay = inv_strray.PrintAlsStrBuilder();//control
            }
            else
            {
                D2 = D1.ToList();
            }
            //=====================================================
            //STEP 2 : Curved Ray Algorithm
            if (B != 0)
            {
                Inv_NetworkWithConstrains inv_network = new Inv_NetworkWithConstrains(_domain, _table, D2, _iP);
                inv_network.SIRT_2();
                _nodesInv = null;
                _nodesInv = inv_network.NodesInv.ToList();
                _sIRT_Res = null;
                _sIRT_Res = inv_network.Sirt_Res.ToList();
            }
            //=====================================================
        }
    }
}
