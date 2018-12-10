using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClassLibrary_TomoGo
{
    public class Input_InversionModel
    {
        //=============================================
        private string _filePath;
        public string FilePath
        {
            get { return _filePath; }
            set { _filePath = value; }
        }
        //=============================================
        private string _fileName;
        public string FileName
        {
            get { return _fileName; }
            set { _fileName = value; }
        }
        //=============================================
        private Inv_SRInfoTable _table;
        public Inv_SRInfoTable Table
        {
            get { return _table; }
            set { _table = value; }
        }
        //=============================================
        public Input_InversionModel copy()
        {
            Input_InversionModel output = new Input_InversionModel();
            output._filePath = string.Copy(this._filePath);
            output._fileName = string.Copy(this._fileName);
            output._table = this._table.copy();
            return output;
        }
    }
}
