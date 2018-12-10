using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

// This code was not orignal written by me (Pengxiang), it´s from somewhere on the internet.
// Just several sentences were modified by me.
// The original comment language is German. I added the corresponding english translation.

namespace ClassLibrary_TomoGo
{
    public class Node
    {
        private string _name;
        public string Name
        {
            get { return _name; }
            set { _name = value; }
        }
        /// <summary>
        /// Konstruktor//Constructor
        /// </summary>
        /// <param name="name">Name des Knotens</param>// the name of the nodes
        /// 
        public Node(string name)
        {
            this._name = name;
        }
        public bool equals(Node n)
        {
            bool eq = false;
            if (_name.Equals(n._name))
                eq = true;
            else
                eq = false;
            return eq;
        }

    }
}
