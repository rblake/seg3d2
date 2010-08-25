/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2009 Scientific Computing and Imaging Institute,
 University of Utah.


 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
 */

#include <Core/Utils/ScopedCounter.h>
#include <Core/Volume/DataVolumeSlice.h>
#include <Core/RenderResources/RenderResources.h>

// Application includes
#include <Application/Tool/ToolFactory.h>
#include <Application/Tools/ThresholdTool.h>
#include <Application/Tools/Actions/ActionThreshold.h>
#include <Application/Tools/detail/MaskShader.h>
#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>
#include <Application/ViewerManager/ViewerManager.h>

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, ThresholdTool )

namespace Seg3D
{

class ThresholdToolPrivate : public Core::Lockable
{
public:
	void handle_threshold_changed();
	void handle_target_layer_changed();
	void handle_seed_points_changed();

	void update_viewers();
	void initialize_gl();

	bool initialized_;
	size_t signal_block_count_;
	ThresholdTool* tool_;
	MaskShaderHandle shader_;
};

void ThresholdToolPrivate::handle_threshold_changed()
{
	if ( this->tool_->target_layer_state_->get() == Tool::NONE_OPTION_C )
	{
		return;
	}
	
	this->update_viewers();
}

void ThresholdToolPrivate::handle_target_layer_changed()
{
	std::string target_layer_id = this->tool_->target_layer_state_->get();
	if ( target_layer_id != Tool::NONE_OPTION_C )
	{
		Core::ScopedCounter signal_block( this->signal_block_count_ );

		DataLayerHandle data_layer = boost::dynamic_pointer_cast< DataLayer >(
			LayerManager::Instance()->get_layer_by_id( target_layer_id ) );
		double min_val = data_layer->get_data_volume()->get_data_block()->get_min();
		double max_val = data_layer->get_data_volume()->get_data_block()->get_max();
		this->tool_->lower_threshold_state_->set_range( min_val, max_val );
		this->tool_->upper_threshold_state_->set_range( min_val, max_val );

		if ( this->tool_->seed_points_state_->get().size() > 0 )
		{
			this->handle_seed_points_changed();
		}
	}
	
	this->update_viewers();
}

void ThresholdToolPrivate::handle_seed_points_changed()
{
	std::string target_layer_id = this->tool_->target_layer_state_->get();
	if ( target_layer_id == Tool::NONE_OPTION_C )
	{
		return;
	}
	
	const Core::StatePointVector::value_type& seed_points = 
		this->tool_->seed_points_state_->get();
	if ( seed_points.size() == 0 )
	{
		this->update_viewers();
		return;
	}
	
	DataLayerHandle data_layer = boost::dynamic_pointer_cast< DataLayer >(
		LayerManager::Instance()->get_layer_by_id( target_layer_id ) );
	Core::DataVolumeHandle data_volume = data_layer->get_data_volume();
	Core::DataBlockHandle data_block = data_volume->get_data_block();
	double min_val = std::numeric_limits< double >::max(); 
	double max_val = std::numeric_limits< double >::min();
	for ( size_t i = 0; i < seed_points.size(); ++i )
	{
		Core::Point pt = data_volume->apply_inverse_grid_transform( seed_points[ i ] );
		int x = Core::Round( pt[ 0 ] );
		int y = Core::Round( pt[ 1 ] );
		int z = Core::Round( pt[ 2 ] );
		if ( x >= 0 && x < static_cast< int >( data_block->get_nx() ) &&
			y >= 0 && y < static_cast< int >( data_block->get_ny() ) &&
			z >= 0 && z < static_cast< int >( data_block->get_nz() ) )
		{
			double val = data_block->get_data_at( static_cast< size_t >( x ),
				static_cast< size_t >( y ), static_cast< size_t >( z ) );
			min_val = Core::Min( min_val, val );
			max_val = Core::Max( max_val, val );
		}
	}
	
	{
		Core::ScopedCounter signal_block( this->signal_block_count_ );
		this->tool_->upper_threshold_state_->set( max_val );
		this->tool_->lower_threshold_state_->set( min_val );
	}

	this->update_viewers();
}

void ThresholdToolPrivate::update_viewers()
{
	if ( this->signal_block_count_ > 0 )
	{
		return;
	}

	ViewerManager::Instance()->update_2d_viewers_overlay();
}

void ThresholdToolPrivate::initialize_gl()
{
	lock_type lock( this->get_mutex() );
	if ( this->initialized_ )
	{
		return;
	}

	{
		Core::RenderResources::lock_type rr_lock( Core::RenderResources::GetMutex() );
		this->shader_.reset( new MaskShader );
		this->shader_->initialize();
	}

	this->shader_->enable();
	this->shader_->set_border_width( 0 );
	this->shader_->set_texture( 0 );
	this->shader_->set_opacity( 0.5f );
	this->shader_->disable();

	this->initialized_ = true;
}

//////////////////////////////////////////////////////////////////////////
// Class Threshold
//////////////////////////////////////////////////////////////////////////

ThresholdTool::ThresholdTool( const std::string& toolid ) :
	SeedPointsTool( Core::VolumeType::DATA_E, toolid ),
	private_( new ThresholdToolPrivate )
{
	this->private_->tool_ = this;
	this->private_->initialized_ = false;
	this->private_->signal_block_count_ = 0;

	this->add_state( "upper_threshold", this->upper_threshold_state_, 5000.0, -5000.0, 5000.0, .01 );
	this->add_state( "lower_threshold", this->lower_threshold_state_, -5000.0, -5000.0, 5000.0, .01 );

	this->private_->handle_target_layer_changed();

	this->add_connection( this->upper_threshold_state_->state_changed_signal_.connect(
		boost::bind( &ThresholdToolPrivate::handle_threshold_changed, this->private_ ) ) );
	this->add_connection( this->lower_threshold_state_->state_changed_signal_.connect(
		boost::bind( &ThresholdToolPrivate::handle_threshold_changed, this->private_ ) ) );
	this->add_connection( this->target_layer_state_->state_changed_signal_.connect(
		boost::bind( &ThresholdToolPrivate::handle_target_layer_changed, this->private_ ) ) );
	this->add_connection( this->seed_points_state_->state_changed_signal_.connect(
		boost::bind( &ThresholdToolPrivate::handle_seed_points_changed, this->private_ ) ) );
}

ThresholdTool::~ThresholdTool()
{
	this->disconnect_all();
}

void ThresholdTool::redraw( size_t viewer_id, const Core::Matrix& proj_mat )
{
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	std::string target_layer_id;
	double min_val, max_val;

	{
		Core::StateEngine::lock_type se_lock( Core::StateEngine::GetMutex() );
		target_layer_id = this->target_layer_state_->get();
		min_val = this->lower_threshold_state_->get();
		max_val = this->upper_threshold_state_->get();
	}

	if ( target_layer_id == Tool::NONE_OPTION_C )
	{
		return;
	}
	
	DataLayerHandle target_layer = boost::dynamic_pointer_cast< DataLayer >(
		LayerManager::Instance()->get_layer_by_id( target_layer_id ) );
	if ( !target_layer )
	{
		CORE_THROW_LOGICERROR( "Target layer '" + target_layer_id + "' doesn't exist" );
	}
	
	Core::DataVolumeSliceHandle target_slice = boost::dynamic_pointer_cast
		< Core::DataVolumeSlice >( viewer->get_volume_slice( target_layer_id ) );
	if ( target_slice->out_of_boundary() )
	{
		return;
	}
	
	this->private_->initialize_gl();

	unsigned int old_tex_unit = Core::Texture::GetActiveTextureUnit();
	Core::Texture::SetActiveTextureUnit( 0 );
	glPushAttrib( GL_TRANSFORM_BIT );
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd( proj_mat.data() );

	std::vector< unsigned char > threshold_mask;
	target_slice->create_threshold_mask( threshold_mask, min_val, max_val, false );
	Core::Texture2DHandle threshold_tex;
	{
		Core::RenderResources::lock_type rr_lock( Core::RenderResources::GetMutex() );
		threshold_tex.reset( new Core::Texture2D );
		threshold_tex->bind();
		threshold_tex->set_min_filter( GL_NEAREST );
		threshold_tex->set_mag_filter( GL_NEAREST );
		threshold_tex->set_image( static_cast< int >( target_slice->nx() ),
			static_cast< int >( target_slice->ny() ), GL_ALPHA, &threshold_mask[ 0 ],
			GL_ALPHA, GL_UNSIGNED_BYTE );
	}

	double voxel_width = ( target_slice->right() - target_slice->left() ) / 
		( target_slice->nx() - 1 );
	double voxel_height = ( target_slice->top() - target_slice->bottom() ) /
		( target_slice->ny() - 1 );
	double left = target_slice->left() - 0.5 * voxel_width;
	double right = target_slice->right() + 0.5 * voxel_width;
	double bottom = target_slice->bottom() - 0.5 * voxel_height;
	double top = target_slice->top() + 0.5 * voxel_height;
	double slice_width = right - left;
	double slice_height = top - bottom;
	Core::Vector slice_x( slice_width, 0.0, 0.0 );
	slice_x = proj_mat * slice_x;
	double slice_screen_width = slice_x.x() / 2.0 * viewer->get_width();
	double slice_screen_height = slice_height / slice_width * slice_screen_width;

	MaskShader::lock_type shader_lock( this->private_->shader_->get_mutex() );
	this->private_->shader_->enable();
	this->private_->shader_->set_pixel_size( static_cast< float >( 1.0 / slice_screen_width ), 
		static_cast< float >( 1.0 /slice_screen_height ) );
	this->private_->shader_->set_color( 1.0f, 1.0f, 0.0f );
	glBegin( GL_QUADS );
		glTexCoord2f( 0.0f, 0.0f );
		glVertex2d( left, bottom );
		glTexCoord2f( 1.0f, 0.0f );
		glVertex2d( right, bottom );
		glTexCoord2f( 1.0f, 1.0f );
		glVertex2d( right, top );
		glTexCoord2f( 0.0f, 1.0f );
		glVertex2d( left, top );
	glEnd();

	threshold_tex->unbind();
	glFinish();
	this->private_->shader_->disable();
	shader_lock.unlock();

	Core::Texture::SetActiveTextureUnit( old_tex_unit );
	glPopMatrix();
	glPopAttrib();

	SeedPointsTool::redraw( viewer_id, proj_mat );
}

void ThresholdTool::execute( Core::ActionContextHandle context )
{
	std::string target_layer;
	double min_val, max_val;
	{
		Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
		target_layer = this->target_layer_state_->get();
		min_val = this->lower_threshold_state_->get();
		max_val = this->upper_threshold_state_->get();
	}

	if ( target_layer != Tool::NONE_OPTION_C )
	{
		ActionThreshold::Dispatch( context, target_layer, min_val, max_val );
	}
}

} // end namespace Seg3D
