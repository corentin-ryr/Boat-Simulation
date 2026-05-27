using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid
{
    public struct ChunkCoord : IEquatable<ChunkCoord>, IComparable<ChunkCoord>
    {
        public int q, r;

        public ChunkCoord(int q, int r) { this.q = q; this.r = r; }

        // A chunk is a flat-top regular hexagon region of the triangle lattice, with
        // circumradius = chunkGridSize · hexRadius (corner lattice point (N,0) lands at
        // world (N·hexRadius, 0)). Chunks tile as a honeycomb at this circumradius, so the
        // axial→world mapping uses the standard flat-top hex layout with size = N·hexRadius.
        // This makes adjacent chunks share an edge and their border lattice vertices coincide.
        public Vector3 WorldCenter(float hexRadius, int chunkGridSize)
        {
            float size = chunkGridSize * hexRadius;

            float x = size * 1.5f * q;
            float z = size * Mathf.Sqrt(3f) * (r + q * 0.5f);

            return new Vector3(x, 0f, z);
        }

        public static ChunkCoord FromWorldPos(Vector3 pos, float hexRadius, int chunkGridSize)
        {
            float size = chunkGridSize * hexRadius;

            float q = 2f / 3f * pos.x / size;

            float r =
                (-1f / 3f * pos.x +
                 Mathf.Sqrt(3f) / 3f * pos.z)
                / size;

            return Round(q, r);
        }

        static ChunkCoord Round(float fq, float fr)
        {
            float fs = -fq - fr;
            int rq = Mathf.RoundToInt(fq), rr = Mathf.RoundToInt(fr), rs = Mathf.RoundToInt(fs);
            float dq = Mathf.Abs(rq - fq), dr = Mathf.Abs(rr - fr), ds = Mathf.Abs(rs - fs);
            if (dq > dr && dq > ds) rq = -rr - rs;
            else if (dr > ds) rr = -rq - rs;
            return new ChunkCoord(rq, rr);
        }

        public IEnumerable<ChunkCoord> HexesInRange(int range)
        {
            for (int dq = -range; dq <= range; dq++)
            {
                int r1 = Mathf.Max(-range, -dq - range);
                int r2 = Mathf.Min(range, -dq + range);
                for (int dr = r1; dr <= r2; dr++)
                    yield return new ChunkCoord(q + dq, r + dr);
            }
        }

        // Lexicographic order on (q, r). Used to pick a single deterministic owner for a
        // dual cell shared by several chunks at a seam, so it is rendered exactly once.
        public int CompareTo(ChunkCoord other) => q != other.q ? q.CompareTo(other.q) : r.CompareTo(other.r);

        public bool Equals(ChunkCoord other) => q == other.q && r == other.r;
        public override bool Equals(object obj) => obj is ChunkCoord c && Equals(c);
        public override int GetHashCode() => q * 73856093 ^ r * 19349663;
        public static bool operator ==(ChunkCoord a, ChunkCoord b) => a.Equals(b);
        public static bool operator !=(ChunkCoord a, ChunkCoord b) => !a.Equals(b);
        public override string ToString() => $"({q},{r})";
    }
}
