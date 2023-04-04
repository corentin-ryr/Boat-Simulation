half3 RGBToHSV(half3 In)
{
    half4 K = half4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    half4 P = lerp(half4(In.bg, K.wz), half4(In.gb, K.xy), step(In.b, In.g));
    half4 Q = lerp(half4(P.xyw, In.r), half4(In.r, P.yzx), step(P.x, In.r));
    half D = Q.x - min(Q.w, Q.y);
    half E = 1e-10;
    return half3(abs(Q.z + (Q.w - Q.y)/(6.0 * D + E)), D / (Q.x + E), Q.x);
}

half3 HSVToRGB(half3 In)
{
    half4 K = half4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    half3 P = abs(frac(In.xxx + K.xyz) * 6.0 - K.www);
    return In.z * lerp(K.xxx, saturate(P - K.xxx), In.y);
}

void HSVLerp_half(half4 A, half4 B, half T, out half4 Out)
{
    A.xyz = RGBToHSV(A.xyz);
    B.xyz = RGBToHSV(B.xyz);

    half t = T; // used to lerp alpha, needs to remain unchanged

    half hue;
    half d = B.x - A.x; // hue difference

    if(A.x > B.x) 
    {
        half temp = B.x;
        B.x = A.x;
        A.x = temp;

        d = -d;
        T = 1-T;
    }

    if(d > 0.5) 
    {
        A.x = A.x + 1;
        hue = (A.x + T * (B.x - A.x)) % 1;
    }

    if(d <= 0.5) hue = A.x + T * d;

    half sat = A.y + T * (B.y - A.y);
    half val = A.z + T * (B.z - A.z);
    half alpha = A.w + t * (B.w - A.w);

    half3 rgb = HSVToRGB(half3(hue,sat,val));
    
    Out = half4(rgb, alpha);
}