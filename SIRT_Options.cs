using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.Statistics;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc: SIRT Algorithm with Options
********************************************/

namespace ClassLibrary_TomoGo
{
    public class SIRT_Options
    {
        public SIRT_Options(Vector<double> b, Matrix<double> A, Vector<double> x, double lamba)
        {
            //Ax=b, lamba is a parameter may used in the algorithm
            _b = b.Clone();
            _A = A.Clone();
            _x = x.Clone();
            _lamda = lamba;
        }
        private Vector<double> _b;
        private Matrix<double> _A;
        private Vector<double> _x;
        private double _lamda;
       //build the matrix M for SIRT-Cimmino
        public Matrix<double> MBuilderCimmino()
        {
            int m = _A.RowCount;
            Vector<double> mDiag = Vector<double>.Build.Dense(m, 0);
            for (int i = 0; i < m; i++)
            {
                double diag = _A.Row(i).L2Norm();
                mDiag[i] = 1 / (m * diag * diag);
            }
            Matrix<double> M = Matrix<double>.Build.DenseOfDiagonalVector(m, m, mDiag);
            return M;
        }
        //build the matrix M for CAV method
        public Matrix<double> MBuilderCAV()
        {
            int m = _A.RowCount;
            int n = _A.ColumnCount;

            // build up Nj
            Vector<double> nDiag = Vector<double>.Build.Dense(n, 0);
            double Nj = 0;
            for (int i = 0; i < n; i++)
            {
                Nj = 0;
                for (int j = 0; j < m; j++)
                {
                    if (_A[j, i] != 0)
                    {
                        Nj = Nj + 1;
                    }
                }
                if (Nj == 0)
                {
                    Nj = 0.01;
                }
                nDiag[i] = Nj;
            }
            // build up M
            Vector<double> mDiag = Vector<double>.Build.Dense(m, 0);
            for (int i = 0; i < m; i++)
            {
                double diag = 0;
                for (int j = 0; j < n; j++)
                {
                    diag = diag + nDiag[j] * _A[i, j] * _A[i, j];
                }
                mDiag[i] = 1 / diag;
            }
            Matrix<double> M = Matrix<double>.Build.DenseOfDiagonalVector(m, m, mDiag);
            return M;
        }
        //build the matrix M for DROP method
        public Matrix<double> MBuilderDROP()
        {
            int m = _A.RowCount;
            Vector<double> mDiag = Vector<double>.Build.Dense(m, 0);
            for (int i = 0; i < m; i++)
            {
                double diag = _A.Row(i).L2Norm();
                mDiag[i] = 1 / (diag * diag);
            }
            Matrix<double> M = Matrix<double>.Build.DenseOfDiagonalVector(m, m, mDiag);
            return M;
        }
        //build the matrix S for DROP method
        public Matrix<double> SBuilderDROP()
        {
            int n = _A.ColumnCount;
            int m = _A.RowCount;
            Vector<double> nDiag = Vector<double>.Build.Dense(n, 1);
            double Nj = 0;
            for (int i = 0; i < n; i++)
            {
                Nj = 0;
                for (int j = 0; j < m; j++)
                {
                    if (_A[j, i] != 0)
                    {
                        //Nj=Nj+1;
                        Nj = Nj + 1;
                    }
                }
                if (Nj == 0)
                    Nj = 0.00001;
                nDiag[i] = 1 / Nj;
            }
            Matrix<double> S = Matrix<double>.Build.DenseOfDiagonalVector(n, n, nDiag);
            return S;
        }
        //the update function for SIT-Cimmino
        public Vector<double> XUpdaterCimmino(Vector<double> XInput, int updateTimes)
        {
            Vector<double> XOutput = XInput.Clone();
            Matrix<double> M = MBuilderCimmino();
            double c = 0;
            for (int i = 0; i < updateTimes; i++)
            {
                Vector<double> r = _b - _A * XInput;
                c = (r * M * r) / (Math.Pow((_A.Transpose() * M * r).L2Norm(), 2.0));
                XOutput = XInput + c * _A.Transpose() * M * r;
                XInput = XOutput.Clone();
            }
            return XOutput;
        }
        //the update function for CAV method with factor
        public Vector<double> XUpdaterCAV(Vector<double> XInput, int updateTimes)
        {
            Vector<double> XOutput = XInput.Clone();
            Matrix<double> M = MBuilderCAV();
            for (int i = 0; i < updateTimes; i++)
            {
                XOutput = XInput + _lamda * _A.Transpose() * M * (_b - _A * XInput);
                XInput = XOutput.Clone();
            }
            return XOutput;
        }
        //the update function for CAV method without factor
        public Vector<double> XUpdaterCAV_NoFactor(Vector<double> XInput, int updateTimes)
        {
            Vector<double> XOutput = XInput.Clone();
            Matrix<double> M = MBuilderCAV();
            for (int i = 0; i < updateTimes; i++)
            {
                XOutput = XInput + _A.Transpose() * M * (_b - _A * XInput);
                XInput = XOutput.Clone();
            }
            return XOutput;
        }
        //the update function for DROP method
        public Vector<double> XUpdaterDROP(Vector<double> XInput, int updateTimes)
        {
            Vector<double> XOutput = XInput.Clone();
            Matrix<double> M = MBuilderDROP();            
            Matrix<double> S = SBuilderDROP();
            for (int i = 0; i < updateTimes; i++)
            {
                Vector<double> r = _b - _A * XInput;
                XOutput = XInput + _lamda * S * _A.Transpose() * M * r;
                XInput = XOutput.Clone();
            }
            return XOutput;
        }
        //the update function for Landweber method
        public Vector<double> XUpdaterLandweber(Vector<double> XInput, int updateTimes)
        {
            double r = 1.0 /Math.Pow( _A.L2Norm(),2.0);
            Vector<double> XOutput = XInput.Clone();
            for (int i = 0; i < updateTimes; i++)
            {
                XOutput = XInput + r * _A.Transpose() * (_b - _A * XInput);
                XInput = XOutput.Clone();
            }
            return XOutput;
        }
        //print function
        public StringBuilder PrintMCimmino() 
        {
            Matrix<double> M = MBuilderCimmino();
            StringBuilder sb = new StringBuilder();
            sb.Append("\r\n------------------------------\r\n");
            sb.Append("M is a diagonal matrix\r\n");
            for (int i = 0; i < M.RowCount; i++)
            {
                sb.Append(M[i, i].ToString("N4") + "   ");
            }
            sb.Append("\r\n------------------------------\r\n");
            return sb;
        }
        //print function 
        public StringBuilder PrintAxb()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("\r\n------------------------------\r\n");
            sb.Append("A in SIRT_Options\r\n");
            for (int i = 0; i < _A.RowCount; i++)
            {
                for (int j = 0; j < _A.ColumnCount; j++)
                {
                    sb.Append(_A[i, j].ToString("N2") + "  ");
                }
                sb.Append("\r\n");
            }
            sb.Append("\r\n------------------------------\r\n");
            sb.Append("x in SIRT_Options\r\n");
            for (int i = 0; i < _x.Count; i++)
                sb.Append(_x[i].ToString("N2") + "    ");
            sb.Append("\r\n------------------------------\r\n");
            sb.Append("b in SIRT_Options\r\n");
            for (int i = 0; i < _b.Count; i++)
                sb.Append(_b[i].ToString("N2") + "    ");
            return sb;
        }


    }
}
