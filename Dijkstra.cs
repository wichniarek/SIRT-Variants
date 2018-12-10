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
    public class Dijkstra
    {
        private double _f;
        public double F
        {
            get { return _f; }
            set { _f = value; }
        }
        /*-------------------------------------------*/
        private List<Node> _nodes;
        public List<Node> Nodes
        {
            get { return _nodes; }
            set { _nodes = value; }
        }
        /*-------------------------------------------*/
        private List<Edge> _edges;
        public List<Edge> Edges
        {
            get { return _edges; }
            set { _edges = value; }
        }
        /*-------------------------------------------*/
        private List<Node> _basis;
        public List<Node> Basis
        {
            get { return _basis; }
            set { _basis = value; }
        }
        /*-------------------------------------------*/
        private Dictionary<string, double> _dist;
        public Dictionary<string, double> Dist
        {
            get { return _dist; }
            set { _dist = value; }
        }
        /*-------------------------------------------*/
        private Dictionary<string, Node> _previous;
        public Dictionary<string, Node> Previous
        {
            get { return _previous; }
            set { _previous = value; }
        }
        /*-------------------------------------------*/
        /// <summary>
        /// Konstruktor//Constructor
        /// </summary>
        /// <param name="edges">Liste aller Kanten</param>//List of all edges
        /// <param name="nodes">Liste aller Knoten</param>//List of all nodes
        /// 
        public Dijkstra(List<Edge> edges, List<Node> nodes)
        {            
            Edges = edges;
            Nodes = nodes;
            Basis = new List<Node>();
            Dist = new Dictionary<string, double>();
            Previous = new Dictionary<string, Node>();

            // Knoten aufnehmen//Add Node 
            foreach (Node n in Nodes)
            {
                Previous.Add(n.Name, null);
                Basis.Add(n);
                Dist.Add(n.Name, double.MaxValue);
            }
        }

        /// <summary>
        /// Knoten zu allen anderen Knoten//Node to the other nodes
        /// </summary>
        /// <param name="start">Startknoten</param>//The begining node
         
        public void calculateDistance(Node start)
        {
            Dist[start.Name] = 0;
            while (Basis.Count > 0)
            {

                Node u = getNodeWithSmallestDistance();

                if (u == null)
                {
                    //MessageBox.Show("here u=null");
                    Basis.Clear();
                }
                else
                {
                    //MessageBox.Show("here u!=null, u="+u.Name);
                    foreach (Node v in getNeighbors(u))
                    {
                        double alt = Dist[u.Name] + getDistanceBetween(u, v);
                        if (alt < Dist[v.Name])
                        {
                            Dist[v.Name] = alt;
                            Previous[v.Name] = u;
                        }
                    }
                    Basis.Remove(u); 
                    

                }
            }
        }
        

        /// <summary>
        /// Liefert den Pfad zum Knoten d//The path to the Node d
        /// </summary>
        /// <param name="d">Zielknote<n/param>//The destination node
        /// <returns></returns>
        public List<Node> getPathTo(Node d)
        {
            List<Node> path = new List<Node>();

            path.Insert(0, d);

            while (Previous[d.Name] != null)
            {
                d = Previous[d.Name];
                path.Insert(0, d);
            }

            return path;
        }

        /// <summary>

        /// </summary>
        /// <returns></returns>
        public Node getNodeWithSmallestDistance()
        {
            double distance = double.MaxValue;
            Node smallest = null;

            foreach (Node n in Basis)
            {
                if (Dist[n.Name] < distance)
                {
                    distance = Dist[n.Name];
                    smallest = n;
                }
            }
            return smallest;
        }

        /// <summary>
        /// Liefert alle Nachbarn von n die noch in der Basis sind//Use all the neigbour nodes of n in the "Basis"
        /// </summary>
        /// <param name="n">Knoten</param>//Nodes
        /// <returns></returns>
        public List<Node> getNeighbors(Node n)
        {
            List<Node> neighbors = new List<Node>();
            foreach (Edge e in Edges)
            {
                if (e.Origin.equals(n))
                {
                    foreach (Node m in Basis)
                    {                        
                        if (m.equals(e.Destination))
                        {
                            neighbors.Add(e.Destination);
                            break;
                        }
                    }
                }

            }
            return neighbors;
        }

        /// <summary>
        /// Liefert die Distanz zwischen zwei Knoten//The distance between two nodes
        /// </summary>
        /// <param name="o">Startknoten</param>//The begining node
        /// <param name="d">Endknoten</param>//The end node
        /// <returns></returns>
        public double getDistanceBetween(Node o, Node d)
        {
            foreach (Edge e in Edges)
            {
                if (e.Origin.equals(o) && e.Destination.equals(d))
                {
                    return e.Distance;
                }
            }
            MessageBox.Show("watch out, function getDistanceBetween() maybe equals 0");
            //what, if nodes o and d are not connected? Answer: then no edge between two nodes:]
            return 0;
        }
    }
}
