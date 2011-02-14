// GLSL vertex shader for slice based volume rendering

uniform bool enable_lighting;
uniform bool enable_fog;
uniform vec3 tex_bbox_min; // Minimum values of texture position in world space
uniform vec3 tex_bbox_size; // Size of texture in world space

varying vec3 light_dir; // Light direction
varying vec3 half_vector; // Half vector between eye and light source. 
                                       // In the case of head light, it's the same as eye vector.
varying float dist; // Distance to light source.

void compute_fog_depth();

void main()
{	
	if ( enable_lighting )
	{
		// The position of the vertex in eye space
		vec4 eye_coord_pos = gl_ModelViewMatrix * gl_Vertex;

		// With headlight, light position is alway at the origin in eye space
		vec3 aux = vec3( -eye_coord_pos );
		light_dir = normalize( aux );
		dist = length(aux);
		half_vector = light_dir;
	}
	
	gl_TexCoord[0].stp = (gl_Vertex.xyz - tex_bbox_min)/tex_bbox_size;
	gl_ClipVertex = gl_ModelViewMatrix * gl_Vertex;
	gl_Position = ftransform();

	if ( enable_fog )
	{
		compute_fog_depth();
	}
} 