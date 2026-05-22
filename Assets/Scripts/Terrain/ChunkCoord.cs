using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid
{
    public struct ChunkCoord : IEquatable<ChunkCoord>
    {
        public int q, r;

        public ChunkCoord(int q, int r) { this.q = q; this.r = r; }

        // Basis vectors for the chunk lattice, derived from the inner hex cell geometry.
        // For a chunk of gridSize=N and hexRadius=R:
        //   e1 (axial step +q): R·√3·(2N+1) in x, 0 in z
        //   e2 (axial step +r): R·√3·N in x,  R·(3N+2) in z
        // These ensure that border vertices of adjacent chunks fall at exactly the
        // same world positions and can be stitched by position deduplication.
        public Vector3 WorldCenter(float hexRadius, int chunkGridSize)
        {
            float e1x = hexRadius * Mathf.Sqrt(3f) * (2 * chunkGridSize + 1);
            float e2x = hexRadius * Mathf.Sqrt(3f) * chunkGridSize;
            float e2z = hexRadius * (3 * chunkGridSize + 2);
            return new Vector3(q * e1x + r * e2x, 0f, r * e2z);
        }

        public static ChunkCoord FromWorldPos(Vector3 pos, float hexRadius, int chunkGridSize)
        {
            float e1x = hexRadius * Mathf.Sqrt(3f) * (2 * chunkGridSize + 1);
            float e2x = hexRadius * Mathf.Sqrt(3f) * chunkGridSize;
            float e2z = hexRadius * (3 * chunkGridSize + 2);
            float fr = pos.z / e2z;
            float fq = (pos.x - fr * e2x) / e1x;
            return Round(fq, fr);
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

        public bool Equals(ChunkCoord other) => q == other.q && r == other.r;
        public override bool Equals(object obj) => obj is ChunkCoord c && Equals(c);
        public override int GetHashCode() => q * 73856093 ^ r * 19349663;
        public static bool operator ==(ChunkCoord a, ChunkCoord b) => a.Equals(b);
        public static bool operator !=(ChunkCoord a, ChunkCoord b) => !a.Equals(b);
        public override string ToString() => $"({q},{r})";
    }
}
