// GLSL vertex shader for per-fragment fog computation

uniform vec2 fog_range;
// Fog depth in eye space. 
// This is a workaround for Macs with ATI X1600/X1900. 
// On those cards accessing gl_FragCoord would cause fallback to software renderer. 
varying float fog_depth;

void compute_fog_depth()
{
	vec4 eye_coord_pos = gl_ModelViewMatrix * gl_Vertex;
	fog_depth = ( -eye_coord_pos.z - fog_range[ 0 ] ) / ( fog_range[ 1 ] - fog_range[ 0 ] );
}
