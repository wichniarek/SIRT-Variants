using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

/********************************************
** auth: Pengxiang Qiu
** date: $time$
** desc: Inversion Parameters
********************************************/

namespace ClassLibrary_TomoGo
{
    public class InversionParameters
    {
        public InversionParameters(int inputNX, int inputNZ)
        {
            NX = inputNX;// the resolution of x direction
            NZ = inputNZ;// the resolution of z direction
            ConBoolMatrix = Matrix<double>.Build.Dense(NZ, NX, 0);// set if this matrix element value is constraint, in other words, if the diffusivity is known, or if the interval if diffusivity is known.
            ConValueMatrix = Matrix<double>.Build.Dense(NZ, NX, InitialD);// this matrix contains the diffusivity distribution
        }
        public int NX;
        public int NZ;
        public int StrRay = 0;//the number of straight ray iteration
        public int CurRay = 0;//the number of curved ray iteration
        public double CriRMS = 0.01;//the stop criteria
        public double InitialD = 1;//the uniform initial D is defined as 1
        public double MinD = double.MinValue;//the lower limit of D
        public double MaxD = double.MaxValue;//the upper limit of D
        public Matrix<double> ConBoolMatrix;
        public Matrix<double> ConValueMatrix;
        public int SIRTOption = 0;//number of SIRTOption indicates the variant of SIRT, for instance, DROP, CAV, Cimmino, etc
        public void UpdateMatrixDim(int inputNX, int inputNZ)
        {
            NX = inputNX;
            NZ = inputNZ;
            ConBoolMatrix = Matrix<double>.Build.Dense(NZ, NX, 0);
            ConValueMatrix = Matrix<double>.Build.Dense(NZ, NX, InitialD);
        }
        public void CopyFrom(InversionParameters IP)
        {
            if (IP.NZ == NZ && IP.NX == NX)
            {
                StrRay = IP.StrRay;
                CurRay = IP.CurRay;
                CriRMS = IP.CriRMS;
                InitialD = IP.InitialD;
                MinD = IP.MinD;
                MaxD = IP.MaxD;
                ConBoolMatrix = IP.ConBoolMatrix;
                ConValueMatrix = IP.ConValueMatrix;
                SIRTOption = IP.SIRTOption;
            }
            else
            {
                MessageBox.Show("Dimension of two inversionparameters are not the same!\r\nError in CopyFrom(InversionParameters)");
            }
        }
        public StringBuilder Show()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("NX="+NX.ToString()+"\n");
            sb.Append("NZ="+NX.ToString()+"\n");
            sb.Append("StrRay="+StrRay.ToString()+"\n");
            sb.Append("CurRay="+CurRay.ToString()+"\n");
            sb.Append("CriRMS="+CriRMS.ToString()+"\n");
            sb.Append("InitialD ="+InitialD .ToString()+"\n");
            sb.Append("MinD="+MinD.ToString()+"\n");
            sb.Append("MaxD="+MaxD.ToString()+"\n");
            sb.Append("ConBoolMatrix:\n"+ConBoolMatrix.ToString()+"\n");
            sb.Append("ConValueMatrix:\n"+ConValueMatrix.ToString()+"\n");
            sb.Append("SIRTOption="+SIRTOption.ToString()+"\n");
            return sb;
        }
    }
}
