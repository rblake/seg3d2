// GLSL vertex shader for rendering a slice
#version 110

uniform bool enable_lighting;
uniform bool enable_fog;

void compute_lighting();
void compute_fog_depth();

void main()
{	
	if ( enable_lighting )
	{
		compute_lighting();
	}
	
	if ( enable_fog )
	{
		compute_fog_depth();
	}
	
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_TexCoord[1] = gl_MultiTexCoord1;
	gl_ClipVertex = gl_ModelViewMatrix * gl_Vertex;
	gl_Position = ftransform();
} 