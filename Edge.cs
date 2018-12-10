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
    public class Edge
    {
        private Node _origin;
        public Node Origin
        {
            get { return _origin; }
            set { _origin = value; }
        }
        /*-------------------------------------------*/
        private Node _destination;
        public Node Destination
        {
            get { return _destination; }
            set { _destination = value; }
        }
        /*-------------------------------------------*/
        private double _distance;
        public double Distance
        {
            get { return _distance; }
            set { _distance = value; }
        }
        /*-------------------------------------------*/
        /// <summary>
        /// Konstruktor//Constructor
        /// </summary>
        /// <param name="origin">Startknoten</param>// the start node
        /// <param name="destination">Zielknoten</param>//the destination node
        /// <param name="distance">Distanz</param>//the distance
        public Edge(Node origin, Node destination, double distance)
        {
            this._origin = origin;
            this._destination = destination;
            this._distance = distance;
        }
        public string print(int n)
        {
            string str = null;
            if (n != 0 && n != 1)
                str = "parameter error\r\nlocation: print() in Edge.cs";
            else
            {
                str = _origin.Name + "  " + _destination.Name + "   " + _distance.ToString() + "\r\n";
                switch (n)
                {
                    case 0:
                        str = "origin" + "  " + "destination" + "   " + "distance" + "\r\n" + str;
                        break;
                    case 1:
                        break;
                }
            }
            return str;
        }

    }
}
