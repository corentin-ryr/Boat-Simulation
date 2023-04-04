// MIT License

// Copyright (c) 2021 NedMakesGames

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Adapted from Ronja Böhringer's tutorials
// https://www.ronja-tutorials.com/post/024-white-noise/
// https://www.ronja-tutorials.com/post/029-tiling-noise/
// Released under CC-BY https://creativecommons.org/licenses/by/2.0/

#ifndef NOISE
#define NOISE

float easeIn(float interpolator) {
	return interpolator * interpolator;
}

float easeOut(float interpolator) {
	return 1 - easeIn(1 - interpolator);
}

float easeInOut(float interpolator) {
	float easeInValue = easeIn(interpolator);
	float easeOutValue = easeOut(interpolator);
	return lerp(easeInValue, easeOutValue, interpolator);
}

float2 modulo(float2 divident, float2 divisor) {
	float2 positiveDivident = divident % divisor + divisor;
	return positiveDivident % divisor;
}

float3 modulo(float3 divident, float3 divisor) {
	float3 positiveDivident = divident % divisor + divisor;
	return positiveDivident % divisor;
}

//to 1d functions

//get a scalar random value from a 3d value
float rand3dTo1d(float3 value, float3 dotDir = float3(12.9898, 78.233, 37.719)) {
	//make value smaller to avoid artefacts
	float3 smallValue = sin(value);
	//get scalar value from 3d vector
	float random = dot(smallValue, dotDir);
	//make value more random by making it bigger and then taking the factional part
	random = frac(sin(random) * 143758.5453);
	return random;
}

float rand2dTo1d(float2 value, float2 dotDir = float2(12.9898, 78.233)) {
	float2 smallValue = sin(value);
	float random = dot(smallValue, dotDir);
	random = frac(sin(random) * 143758.5453);
	return random;
}

float rand1dTo1d(float3 value, float mutator = 0.546) {
	float random = frac(sin(value + mutator) * 143758.5453);
	return random;
}

//to 2d functions

float2 rand3dTo2d(float3 value) {
	return float2(
		rand3dTo1d(value, float3(12.989, 78.233, 37.719)),
		rand3dTo1d(value, float3(39.346, 11.135, 83.155))
	);
}

float2 rand2dTo2d(float2 value) {
	return float2(
		rand2dTo1d(value, float2(12.989, 78.233)),
		rand2dTo1d(value, float2(39.346, 11.135))
	);
}

float2 rand1dTo2d(float value) {
	return float2(
		rand2dTo1d(value, 3.9812),
		rand2dTo1d(value, 7.1536)
	);
}

//to 3d functions

float3 rand3dTo3d(float3 value) {
	return float3(
		rand3dTo1d(value, float3(12.989, 78.233, 37.719)),
		rand3dTo1d(value, float3(39.346, 11.135, 83.155)),
		rand3dTo1d(value, float3(73.156, 52.235, 09.151))
	);
}

float3 rand2dTo3d(float2 value) {
	return float3(
		rand2dTo1d(value, float2(12.989, 78.233)),
		rand2dTo1d(value, float2(39.346, 11.135)),
		rand2dTo1d(value, float2(73.156, 52.235))
	);
}

float3 rand1dTo3d(float value) {
	return float3(
		rand1dTo1d(value, 3.9812),
		rand1dTo1d(value, 7.1536),
		rand1dTo1d(value, 5.7241)
	);
}

// Shader graph custom functions
void Rand3dTo1d_float(float3 Value, out float Out) {
	rand3dTo1d(Value, Out);
}

void Rand2dTo1d_float(float2 Value, out float Out) {
	rand2dTo1d(Value, Out);
}

void Rand1dTo1d_float(float3 Value, out float Out) {
	rand1dTo1d(Value, Out);
}

void Rand3dTo2d_float(float3 value, out float2 Out) {
	Out = rand3dTo2d(value);
}

void Rand2dTo2d_float(float2 value, out float2 Out) {
	Out = rand2dTo2d(value);
}

void Rand1dTo2d_float(float value, out float2 Out) {
	Out = rand1dTo2d(value);
}

void Rand3dTo3d_float(float3 value, out float3 Out) {
	Out = rand3dTo3d(value);
}

void Rand2dTo3d_float(float2 value, out float3 Out) {
	Out = rand2dTo3d(value);
}

void Rand1dTo3d_float(float value, out float3 Out) {
	Out = rand1dTo3d(value);
}

// Perlin / gradient noise

float perlinNoise(float2 value, float2 period) {
	float2 cellsMimimum = floor(value);
	float2 cellsMaximum = ceil(value);

	cellsMimimum = modulo(cellsMimimum, period);
	cellsMaximum = modulo(cellsMaximum, period);

				//generate random directions
	float2 lowerLeftDirection = rand2dTo2d(float2(cellsMimimum.x, cellsMimimum.y)) * 2 - 1;
	float2 lowerRightDirection = rand2dTo2d(float2(cellsMaximum.x, cellsMimimum.y)) * 2 - 1;
	float2 upperLeftDirection = rand2dTo2d(float2(cellsMimimum.x, cellsMaximum.y)) * 2 - 1;
	float2 upperRightDirection = rand2dTo2d(float2(cellsMaximum.x, cellsMaximum.y)) * 2 - 1;

	float2 fraction = frac(value);

				//get values of cells based on fraction and cell directions
	float lowerLeftFunctionValue = dot(lowerLeftDirection, fraction - float2(0, 0));
	float lowerRightFunctionValue = dot(lowerRightDirection, fraction - float2(1, 0));
	float upperLeftFunctionValue = dot(upperLeftDirection, fraction - float2(0, 1));
	float upperRightFunctionValue = dot(upperRightDirection, fraction - float2(1, 1));

	float interpolatorX = easeInOut(fraction.x);
	float interpolatorY = easeInOut(fraction.y);

				//interpolate between values
	float lowerCells = lerp(lowerLeftFunctionValue, lowerRightFunctionValue, interpolatorX);
	float upperCells = lerp(upperLeftFunctionValue, upperRightFunctionValue, interpolatorX);

	float noise = lerp(lowerCells, upperCells, interpolatorY);
	return noise;
}

float sampleLayeredNoise(float2 value, float period, float persistance, float roughness, float octaves) {
	float noise = 0;
	float frequency = 1;
	float factor = 1;

	for(int i = 0; i < octaves; i++) {
		noise = noise + perlinNoise(value * frequency + i * 0.72354, period * frequency) * factor;
		factor *= persistance;
		frequency *= roughness;
	}

	return noise;
}

void PerlinNoise_float(float2 Value, float2 Period, out float Out) {
	Out = perlinNoise(Value, Period) + 0.5;
}

void PerlinNoise_float(float2 Value, float2 Period, float Persistance, float Roughness, float Octaves, out float Out) {
	Out = sampleLayeredNoise(Value, Period, Persistance, Roughness, Octaves) + 0.5;

}

// Voronoi noise

float3 voronoiNoise(float3 value, float3 period) {
	float3 baseCell = floor(value);

	//first pass to find the closest cell
	float minDistToCell = 10;
	float3 toClosestCell;
	float3 closestCell;
	[unroll]
	for(int x1 = -1; x1 <= 1; x1++) {
		[unroll]
		for(int y1 = -1; y1 <= 1; y1++) {
			[unroll]
			for(int z1 = -1; z1 <= 1; z1++) {
				float3 cell = baseCell + float3(x1, y1, z1);
				float3 tiledCell = modulo(cell, period);
				float3 cellPosition = cell + rand3dTo3d(tiledCell);
				float3 toCell = cellPosition - value;
				float distToCell = length(toCell);
				if(distToCell < minDistToCell) {
					minDistToCell = distToCell;
					closestCell = cell;
					toClosestCell = toCell;
				}
			}
		}
	}

	//second pass to find the distance to the closest edge
	float minEdgeDistance = 10;
	[unroll]
	for(int x2 = -1; x2 <= 1; x2++) {
		[unroll]
		for(int y2 = -1; y2 <= 1; y2++) {
			[unroll]
			for(int z2 = -1; z2 <= 1; z2++) {
				float3 cell = baseCell + float3(x2, y2, z2);
				float3 tiledCell = modulo(cell, period);
				float3 cellPosition = cell + rand3dTo3d(tiledCell);
				float3 toCell = cellPosition - value;

				float3 diffToClosestCell = abs(closestCell - cell);
				bool isClosestCell = diffToClosestCell.x + diffToClosestCell.y + diffToClosestCell.z < 0.1;
				if(!isClosestCell) {
					float3 toCenter = (toClosestCell + toCell) * 0.5;
					float3 cellDifference = normalize(toCell - toClosestCell);
					float edgeDistance = dot(toCenter, cellDifference);
					minEdgeDistance = min(minEdgeDistance, edgeDistance);
				}
			}
		}
	}

	float random = rand3dTo1d(closestCell);
	return float3(minDistToCell, random, minEdgeDistance);
}

void VoronoiNoise_float(float3 Value, float3 Period, out float MinCellDistance, out float Random, out float MinEdgeDistance) {
	float3 ret = voronoiNoise(Value, Period);
	MinCellDistance = ret.x;
	Random = ret.y;
	MinEdgeDistance = ret.z;
}

#endif