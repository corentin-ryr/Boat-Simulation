using System;
using UnityEngine;

namespace Geometry
{
    public struct Edge
    {
        public int X;
        public int Y;

        public Edge(int x, int y)
        {
            this.X = x;
            this.Y = y;
        }

        public override bool Equals(object obj)
        {
            if (obj is Edge) return this.Equals((Edge)obj);
            else return false;
        }

        public bool Equals(Edge other)
        {
            return (((this.X == other.X) && (this.Y == other.Y)) || ((this.X == other.Y) && (this.Y == other.X)));
        }

        public override int GetHashCode()
        {
            return (Math.Min(this.X, this.Y).GetHashCode() + Math.Max(this.X, this.Y).GetHashCode());
        }

        public override string ToString()
        {
            return string.Format("{{X:{0} Z:{1}}}", this.X, this.Y);
        }
    }


}