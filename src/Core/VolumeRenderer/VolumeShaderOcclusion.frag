// GLSL fragment shader for rendering a slice

uniform sampler3D vol_tex;
uniform sampler1D diffuse_lut;
uniform sampler1D occlusion_sample_lut;
uniform sampler2D occlusion_buffer;

// Number of samples contained in the occlusion_sample_lut
uniform int num_of_occlusion_samples;
// Occlusion extent (radius of sampling disk) in world space
uniform float occlusion_extent;
// Slice distance normalized by volume unit length
uniform float normalized_slice_distance; 

varying vec2 correction_factor;
varying vec4 clip_space_pos;

float volume_lookup( vec3 tex_coord )
{
	float val = texture3D( vol_tex, tex_coord ).a;
	return val;
}

void main()
{
	float voxel_val = volume_lookup( gl_TexCoord[0].stp );
	vec4 diffuse_color = texture1D( diffuse_lut, voxel_val );

	float transparency = pow( 1.0 - diffuse_color.a, normalized_slice_distance );
	//float transparency = exp( -diffuse_color.a * normalized_slice_distance );
	float alpha = 1.0 - transparency;

	// Render to the eye buffer
	vec2 device_space_pos = clip_space_pos.xy / clip_space_pos.w;
	vec2 tex_space_pos = device_space_pos * 0.5 + vec2( 0.5 );
	float occlusion = texture2D( occlusion_buffer, tex_space_pos ).r;
	gl_FragData[ 0 ] = vec4( diffuse_color.rgb * occlusion * alpha, alpha );

	// Render to the next occlusion buffer
	occlusion = 0.0;
	for ( int i = 0; i < num_of_occlusion_samples; ++i )
	{
		vec2 sample_pos = texture1D( occlusion_sample_lut, 
			( float( i ) + 0.5 ) / float( num_of_occlusion_samples ) ).xy;
		sample_pos *= occlusion_extent;
		sample_pos = device_space_pos + correction_factor * sample_pos;
		sample_pos = sample_pos * 0.5 + vec2( 0.5 );
		occlusion += texture2D( occlusion_buffer, sample_pos ).r;
	}
	occlusion = occlusion / float( num_of_occlusion_samples ) * transparency;
	gl_FragData[ 1 ] = vec4( occlusion );
}
