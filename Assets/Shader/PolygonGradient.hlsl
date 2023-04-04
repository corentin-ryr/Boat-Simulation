// Licensed under the Unity Companion Package License v1.0
// https://docs.unity3d.com/Packages/com.unity.shadergraph@7.1/license/LICENSE.html

#ifndef POLYGON_GRADIENT
#define POLYGON_GRADIENT

void PolygonGradient_float(float2 UV, float Sides, float Width, float Height, out float Out) {
	float pi = PI;
	float aWidth = Width * cos(pi / Sides);
	float aHeight = Height * cos(pi / Sides);
	float2 uv = (UV * 2 - 1) / float2(aWidth, aHeight);
	uv.y *= -1;
	float pCoord = atan2(uv.x, uv.y);
	float r = 2 * pi / Sides;
	float distance = cos(floor(0.5 + pCoord / r) * r - pCoord) * length(uv);
	Out = 1 - distance;
}

#endif