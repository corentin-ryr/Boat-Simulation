using System;
using System.Collections.Generic;
using UnityEngine;

namespace TerrainGrid.MeshBuilders
{
    // Lighthouse = stacked silhouette inspired by Petit Minou (Brittany):
    //
    //         ▲      finial mast
    //        ╱ ╲     dome (cone, lantern-room cap)
    //       │   │    lantern room (slightly narrower cylinder)
    //      ┌─────┐   gallery ring (wider band — the balcony)
    //      │     │
    //      │     │   stone shaft (gently tapered cylinder, the bulk)
    //      │     │
    //      └─────┘
    //       ▔▔▔▔▔    harbour deck (rendered by DockMesh at Y = PlatformTopY)
    //
    // Sized to read as a real landmark at game scale: the shaft fills
    // most of the dual cell footprint and the whole tower reaches ~15 m
    // above the deck (≈5× a Building), so it's visible from across the
    // bay. Single material today — the iconic red lantern + dome are
    // shape-only for now; introducing a second material would mean
    // splitting the tower into a second overlay (see "Out of scope" note
    // at the bottom of this file).
    //
    // Anchored on the shared harbour deck — Y0 = DockMesh.PlatformTopY.
    public static class LighthouseMesh
    {
        // Must match DockMesh.PlatformTopY — the shaft sits on the deck
        // that DockMesh draws for this same cell.
        const float PlatformTopY = 2f;

        // Silhouette layers (radii = circumradius / centre-to-corner in
        // world units; height in metres). With Sides = 6 the visible
        // edge-midpoint distance (apothem) is r·√3/2 ≈ 0.866·r.
        const float ShaftBaseRadius = 1.10f;   // pushes near the dual-cell corners
        const float ShaftTopRadius  = 1.00f;
        const float ShaftHeight     = 9.0f;

        const float GalleryRadius   = 1.20f;   // overhangs the shaft → balcony band
        const float GalleryHeight   = 0.5f;

        const float LanternRadius   = 0.90f;   // the glass-walled lamp room
        const float LanternHeight   = 2.0f;

        const float DomeBaseRadius  = 0.92f;
        const float DomePeakRadius  = 0.12f;
        const float DomeHeight      = 1.8f;

        const float FinialBottomR   = 0.06f;
        const float FinialTopR      = 0.0f;
        const float FinialHeight    = 1.0f;

        // 6 sides for a hexagonal silhouette — matches the underlying
        // primal hex lattice and reads as masonry rather than a turned
        // cylinder.
        const int Sides = 6;

        public static Mesh Build(PrimalChunk primal, Func<ChunkCoord, PrimalChunk> lookupPrimal)
        {
            _ = lookupPrimal;

            var verts = new List<Vector3>();
            var tris = new List<int>();

            TileProperty[] tp = primal.TileProperties;
            Vector3[] centroids = primal.PrimalCellCentroids;
            if (tp == null || centroids == null) return MeshPrimitives.EmptyMesh();

            for (int i = 0; i < tp.Length && i < centroids.Length; i++)
            {
                if (tp[i].Kind != TileKind.Lighthouse) continue;
                Vector3 c = centroids[i];

                // Stack Y-coordinates. All anchored on the fixed deck,
                // never on the cell centroid, so every Lighthouse cell
                // has identical silhouette regardless of seabed depth.
                float y0 = PlatformTopY;
                float y1 = y0 + ShaftHeight;
                float y2 = y1 + GalleryHeight;
                float y3 = y2 + LanternHeight;
                float y4 = y3 + DomeHeight;
                float y5 = y4 + FinialHeight;

                // Stone shaft — slight inward taper top→bottom for stability read.
                MeshPrimitives.AddTaperedPrism(verts, tris, c.x, c.z,
                    y0, ShaftBaseRadius, y1, ShaftTopRadius, Sides);

                // Gallery ring — wider than shaft, short height; the
                // band that visually anchors the lantern room. Built as
                // a flat-walled prism (rB == rT) so it reads as a solid
                // balcony rather than a taper.
                MeshPrimitives.AddTaperedPrism(verts, tris, c.x, c.z,
                    y1, GalleryRadius, y2, GalleryRadius, Sides);

                // Lantern room — narrower than the gallery → visible
                // setback above the balcony.
                MeshPrimitives.AddTaperedPrism(verts, tris, c.x, c.z,
                    y2, LanternRadius, y3, LanternRadius, Sides);

                // Dome — broad base tapering nearly to a point; the
                // iconic curved cap. Approximated as a single cone; a
                // proper hemispherical dome would need a stack of rings
                // and isn't worth the extra triangles at this distance.
                MeshPrimitives.AddTaperedPrism(verts, tris, c.x, c.z,
                    y3, DomeBaseRadius, y4, DomePeakRadius, Sides);

                // Finial mast — thin spike on top of the dome. Tapers to
                // a point for a clean silhouette.
                MeshPrimitives.AddTaperedPrism(verts, tris, c.x, c.z,
                    y4, FinialBottomR, y5, FinialTopR, Sides);
            }

            return MeshPrimitives.Finalize(verts, tris);
        }

        // Out of scope (for now):
        //   • Red lantern room + red dome cap. The lighthouse overlay
        //     has one material (LighthouseMaterial). Two-tone would
        //     need either a second overlay (top-half draws with a
        //     separate red material) or per-vertex colour into a tinted
        //     URP/Lit material. Both are a half-day refactor — done
        //     when the silhouette is locked in.
        //   • Annex building at the base. The reference has a low keeper's
        //     cottage attached; would attach to one wedge boundary of
        //     the cell. Deferred until we know if the cell's orientation
        //     ("which side faces water") is exposed.
        //   • Window/light emission. Needs a lit shader pass — out of
        //     scope for the URP/Lit material baseline.
    }
}
