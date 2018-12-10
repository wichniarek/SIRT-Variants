using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Statistics;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc:The Straight Ray Inversion with Constrains
********************************************/

namespace ClassLibrary_TomoGo
{
    public class Inv_StraightRayWithConstrain
    {
        public Inv_StraightRayWithConstrain(float[] domain, Inv_SRInfoTable table, List<double> d, InversionParameters inputIP)
        {

            if (domain.Length != 4 || inputIP.NX < 1 || inputIP.NZ < 1)
            {
                MessageBox.Show("input parameter error in Inv_StraightRay!");
                //check if the input
            }
            _domain = new float[4];
            _domain[0] = domain[0]; _domain[1] = domain[1]; _domain[2] = domain[2]; _domain[3] = domain[3];

            _table = table.copy();
            _iP = new InversionParameters(inputIP.NX, inputIP.NZ);
            _iP.CopyFrom(inputIP);
            conMinX = Vector<double>.Build.Dense(_iP.NX * _iP.NZ, double.MinValue);
            conMaxX = Vector<double>.Build.Dense(_iP.NX * _iP.NZ, double.MaxValue);

            CreateConValues(_iP);// build up conMinValue and conMaxValue
            //MessageBox.Show("After\r\n" + conMaxX.ToString() + "\r\n" + conMinX.ToString());
            _initialDiffusivity = d.ToList();// these initial D in cells are mightbe different
            _tMeasured = TsGenerator1(_table);// initial travel time
            _tProcessed = TsGenerator2(_table);// travel time after transition
            _fRs = DsGenerator(_domain, _iP.NX, _iP.NZ);// construt the elements (rectangles)
            _result = new List<SIRT_Result_StraightRay>();
        }
        private InversionParameters _iP;
        private Vector<double> conMinX;
        private Vector<double> conMaxX;
        private float[] _domain;
        private Inv_SRInfoTable _table;
        private List<double> _tMeasured;
        private List<double> _tProcessed;
        private List<double> _initialDiffusivity;
        private List<FRechtangle> _fRs;
        private List<SIRT_Result_StraightRay> _result;
        public List<SIRT_Result_StraightRay> Result
        {
            get { return _result; }
            set { _result = value; }
        }
        //=============================================
        private List<double> TsGenerator1(Inv_SRInfoTable table)
        {
            List<double> output = new List<double>();
            foreach (Inv_SRInfoLine L in table.SRInfoLineList)
            {
                Inv_SRInfoLine ttl = L.copy();
                ttl.TravelTime = ttl.TravelTime;
                output.Add(ttl.TravelTime);
            }
            return output;
        }
        private List<double> TsGenerator2(Inv_SRInfoTable table)
        {
            List<double> output = new List<double>();
            foreach (Inv_SRInfoLine L in table.SRInfoLineList)
            {
                Inv_SRInfoLine ttl = L.copy();
                ttl.TravelTime = (float)Math.Pow(6 * table.F_alpha_d * ttl.TravelTime, 0.5);
                output.Add(ttl.TravelTime);
            }
            return output;
        }
        public List<FRechtangle> DsGenerator(float[] domain, int nx, int nz)
        {
            List<FRechtangle> List = new List<FRechtangle>();
            int elementcount = nx * nz;
            float elementlength = (domain[1] - domain[0]) / nx;
            float elementheight = (domain[3] - domain[2]) / nz;
            for (int i = 0; i < elementcount; i++)
            {
                int row = i / nx;//show the row where D(i) in;
                int col = i % nx;//show the column where D(i) in;                
                FRechtangle fr = new FRechtangle(domain[0] + col * elementlength, domain[0] + (col + 1) * elementlength,
                                                 domain[3] - (row + 1) * elementheight, domain[3] - row * elementheight);
                fr.Name = i.ToString();
                List.Add(fr);
            }
            return List;
        }
        public void StRayAlgo()//straight ray algorithm
        {
            if (_result.Count == 0)
            {
                _result.Clear();
            }
            int n_d = _initialDiffusivity.Count;//he number of cells, it equals the number of columns of matrix A
            int n_t = _tProcessed.Count;//the number of travel times, it equals the number of rows of matrix A
            // b=Ax, x->x2   
            // consider b as time, A as distance, x as speed
            double xMax = 1.0 / (Math.Sqrt(_iP.MinD));// set the constraint, the max of x, x=1/(root of diffusivity)
            double xMin = 1.0 / (Math.Sqrt(_iP.MaxD));// set the constraint, the min of x, x=1/(root of diffusivity)
            List<double> estimatedD = _initialDiffusivity.ToList();
            Vector<double> x = Vector<double>.Build.Dense(n_d, 1);
            for (int i = 0; i < n_d; i++)
            {
                x[i] = 1.0 / (Math.Sqrt(estimatedD[i]));//the speed vector, x=1/(root of diffusivity)
            }

            #region
            int iter = 0;
            double RMS = _iP.CriRMS + 1;
            while (iter < _iP.StrRay & RMS > _iP.CriRMS)
            {
                //STEP 1 
                SIRT_Result_StraightRay Record = new SIRT_Result_StraightRay();
                Record.Iter = iter;
                Record.OldProcessedD = x.ToList();
                //====================================================
                //STEP 2 : Set matrix A
                Matrix<double> A = Matrix<double>.Build.Dense(n_t, n_d, 0);
                for (int row = 0; row < n_t; row++)
                {
                    Inv_SRInfoLine infoLine = (Inv_SRInfoLine)_table.SRInfoLineList[row];
                    MathNet.Spatial.Euclidean.Point2D PStart = new MathNet.Spatial.Euclidean.Point2D(infoLine.Start.Coor[0], infoLine.Start.Coor[2]);
                    MathNet.Spatial.Euclidean.Point2D PEnd = new MathNet.Spatial.Euclidean.Point2D(infoLine.End.Coor[0], infoLine.End.Coor[2]);
                    MathNet.Spatial.Euclidean.Line2D Traj = new MathNet.Spatial.Euclidean.Line2D(PStart, PEnd);
                    for (int col = 0; col < n_d; col++)
                    {
                        CohenSutherland coh = new CohenSutherland();
                        List<MathNet.Spatial.Euclidean.Point2D> Ps = coh.CohenSutherlandLineClip(_fRs[col], PStart, PEnd).ToList();
                        if (Ps.Count > 1)
                        {
                            A[row, col] = Ps[0].DistanceTo(Ps[1]);
                        }

                    }
                }
                Record.A = A.Clone();
                Record.ProcessedTime = (A * x).ToList();
                //====================================================
                //STEP 3 : set the start node and the end node from Record
                foreach (object obj in _table.SRInfoLineList)
                {
                    Inv_SRInfoLine infoLine = (Inv_SRInfoLine)obj;
                    Record.Start.Add(infoLine.Start);
                    Record.End.Add(infoLine.End);
                }
                //====================================================
                //STEP 4 : cimmino calculation    
                Vector<double> b = Vector<double>.Build.DenseOfArray(_tProcessed.ToArray());
                SIRT_Options so = new SIRT_Options(b, A, x, 1);// b=Ax
                Vector<double> x2 = so.XUpdaterDROP(x, 1);
                //====================================================
                //STEP 7 : update x and diffusivity for next iteration
                x = x2.Clone();
                for (int i = 0; i < n_d; i++)
                {
                    if (x[i] < conMinX[i])
                        x[i] = conMinX[i];
                    if (x[i] > conMaxX[i])
                        x[i] = conMaxX[i];
                }
                Record.NewProcessedD = x.ToList();
                //====================================================
                _result.Add(Record);//output
                iter = iter + 1;
            }
            #endregion

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
        public System.Text.StringBuilder PrintAlsStrBuilder()
        {
            System.Text.StringBuilder sbPrint = new System.Text.StringBuilder();
            if (Result != null)
            {
                foreach (SIRT_Result_StraightRay R in Result)
                {
                    sbPrint.Append(R.ResultAlsStrBuilder());
                }
            }
            return sbPrint;
        }
    }
}
