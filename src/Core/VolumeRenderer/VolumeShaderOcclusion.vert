// GLSL vertex shader for slice based volume rendering

uniform vec3 tex_bbox_min; // Minimum values of texture position in world space
uniform vec3 tex_bbox_size; // Size of texture in world space

varying vec2 correction_factor;
varying vec4 clip_space_pos;

void main()
{
	// The position of the vertex in eye space
	vec4 eye_coord_pos = gl_ModelViewMatrix * gl_Vertex;
	vec4 proj_scale = gl_ProjectionMatrix * vec4( 1.0, 1.0, eye_coord_pos.z, 1.0 );
	correction_factor = vec2( proj_scale.x / proj_scale.w, proj_scale.y / proj_scale.w );

	gl_TexCoord[0].stp = (gl_Vertex.xyz - tex_bbox_min)/tex_bbox_size;
	gl_ClipVertex = gl_ModelViewMatrix * gl_Vertex;
	clip_space_pos = ftransform();
	gl_Position = clip_space_pos;
} 