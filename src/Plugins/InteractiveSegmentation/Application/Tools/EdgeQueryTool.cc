/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2013 Scientific Computing and Imaging Institute,
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

// Core includes
#include <Core/Volume/VolumeSlice.h>
#include <Core/Volume/MaskVolumeSlice.h>

#include <Core/State/Actions/ActionSetAt.h>
#include <Core/State/Actions/ActionSet.h>
#include <Core/Interface/Interface.h>

#include <Core/DataBlock/SliceType.h>

#include <Core/Utils/Lockable.h>

// Application includes
#include <Application/Tool/ToolFactory.h>
#include <Application/Layer/Layer.h>
#include <Application/Layer/LayerManager.h>
#include <Application/Viewer/Viewer.h>
#include <Application/ViewerManager/ViewerManager.h>

#include <Plugins/InteractiveSegmentation/Application/Tools/Actions/ActionEdgeQuery.h>

#include <Plugins/InteractiveSegmentation/Application/Tools/EdgeQueryTool.h>

// test
#include <iostream>
// test

// Register the tool into the tool factory
SCI_REGISTER_TOOL( Seg3D, EdgeQueryTool )

namespace Seg3D
{

  
//////////////////////////////////////////////////////////////////////////
// Class EdgeQuery
//////////////////////////////////////////////////////////////////////////

#define EDGE_QUERY_SIZE 2  

struct LineSegment
{
  LineSegment(const Core::Point& p1, const Core::Point& p2, const int label)
    : p1_(p1), p2_(p2), label_(label), selected_(true), hovering_(false), notInSlice_(false)
  {}  
  
  Core::Point p1_;
  Core::Point p2_;
  
  int label_;
  
  bool selected_;
  bool hovering_;
  bool notInSlice_;
};

class EdgeQuery
{
public:
  EdgeQuery() {}

  // TODO: implement...
  bool is_valid() const
  {
    return lines_.size() >= EDGE_QUERY_SIZE;
  }

  //  bool invalidate()  { return false; }
  
  void create_edgeQuery(const Core::Point& p1, const Core::Point& p2, const Core::Point& p3)
  {
    lines_.reserve(EDGE_QUERY_SIZE);
    
    LineSegment e1(p1, p2, 0);
    LineSegment e2(p2, p3, 1);
    lines_.push_back(e1);
    lines_.push_back(e2);
  }
  
  bool allUnselected() const
  {
    return (! lines_[0].selected_) && (! lines_[1].selected_);
  }
  
  bool allSelected() const
  {
    return ( lines_[0].selected_ && lines_[1].selected_ );
  }
  
  bool selected(int index)
  {
    if (! validIndex(index) )
    {
      return false;
    }
    return lines_[index].selected_;
  }
  
  void select(int index)
  {
    if (! validIndex(index) )
    {
      return;
    }
    lines_[index].selected_ = true;
  }
  
  void unselect(int index)
  {
    if (! validIndex(index) )
    {
      return;
    }
    lines_[index].selected_ = false;
  }

  void hover(int index)
  {
    if (! validIndex(index) )
    {
      return;
    }
    lines_[index].hovering_ = true;
  }
  
  void unsetHover(int index)
  {
    if (! validIndex(index) )
    {
      return;
    }
    lines_[index].hovering_ = false;
  }
  
  bool hovering(int index)
  {
    if (! validIndex(index) )
    {
      return false;
    }
    return lines_[index].hovering_;
  }

  void clearHover()
  {
    for (int i = 0; i < EDGE_QUERY_SIZE; ++i)
    {
      lines_[i].hovering_ = false;
    }
  }
  
  LineSegment getEdge(int index) const
  {
    // TODO: should have error checking here...exception?
    // this implementation assumes validity already checked

    return lines_[index];
  }
  
private:
  bool validIndex(int index)
  {
    if (index > EDGE_QUERY_SIZE - 1)
    {
      // TODO: should throw exception?
      return false;
    }
    return true;
  }

  std::vector<LineSegment> lines_;
};


//////////////////////////////////////////////////////////////////////////
// Class EdgeQueryToolPrivate
//////////////////////////////////////////////////////////////////////////

class EdgeQueryToolPrivate
{
public:
	void handle_vertices_changed();
  void handle_selection_changed();
//  void handle_save_changed();

	void execute( Core::ActionContextHandle context, bool erase, bool save,
		ViewerHandle viewer = ViewerHandle() );
  bool point_in_slice( ViewerHandle viewer, const Core::Point& world_point, 
    double& point_depth ) const;

	EdgeQueryTool* tool_;
  EdgeQuery edgeQuery_;
};

bool EdgeQueryToolPrivate::point_in_slice( ViewerHandle viewer, const Core::Point& world_point, 
                                            double& point_depth ) const
{
  // May be called from interface or rendering thread
  
  // Basically in_slice has to be within epsilon of this slice so that editing won't move the
  // measurement points
  Core::VolumeSliceHandle volume_slice = viewer->get_active_volume_slice();
  if ( ! volume_slice )
  {
    return false;
  }
  
  double slice_depth = volume_slice->depth();
  double i_pos, j_pos;
  volume_slice->project_onto_slice( world_point, i_pos, j_pos, point_depth );

  //std::cerr << "point_depth=" << point_depth << ", slice_depth=" << slice_depth << std::endl;

  if ( point_depth == slice_depth )
  {
    return true;
  }
  return false;
}
  
void EdgeQueryToolPrivate::handle_vertices_changed()
{
  //Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );

  std::vector< Core::Point > points = this->tool_->vertices_state_->get();
  if (points.size() >= 3)
  {
    this->edgeQuery_.create_edgeQuery(points[0], points[1], points[2]);
 
    Core::ActionSet::Dispatch( Core::Interface::GetMouseActionContext(),
                              this->tool_->vertices_state_, points );

    for (int i = 0; i < EDGE_QUERY_SIZE; ++i)
    {
      if ( this->edgeQuery_.selected(i) )
      {
        this->tool_->selectedEdges_state_->add(1);
      }
      else
      {
        this->tool_->selectedEdges_state_->add(0);
      }
    }

    ViewerManager::Instance()->update_2d_viewers_overlay();
    
//    {
//      Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
//      this->save_state_->set(false);
//    }

  }
}

void EdgeQueryToolPrivate::handle_selection_changed()
{
	ViewerManager::Instance()->update_2d_viewers_overlay();
}
	
void EdgeQueryToolPrivate::execute( Core::ActionContextHandle context, 
								  bool save, bool stop, ViewerHandle viewer )
{
	Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );

	if ( !this->tool_->valid_target_state_->get() )
	{
		return;
	}
	
	// If no viewer specified, use the current active viewer
	if ( !viewer )
	{
		int active_viewer = ViewerManager::Instance()->active_viewer_state_->get();
		if ( active_viewer < 0 )
		{
			return;
		}
		viewer = ViewerManager::Instance()->get_viewer( static_cast< size_t >( active_viewer ) );
	}
	
	if ( !viewer || viewer->is_volume_view() )
	{
		return;
	}
	
	Core::MaskVolumeSliceHandle volume_slice = boost::dynamic_pointer_cast
		< Core::MaskVolumeSlice >( viewer->get_volume_slice( 
		this->tool_->target_layer_state_->get() ) );
	if ( !volume_slice )
	{
		return;
	}

	LayerHandle layer = LayerManager::Instance()->find_layer_by_id( 
		this->tool_->target_layer_state_->get() );
	if ( !layer->is_visible( viewer->get_viewer_id() ) || layer->locked_state_->get() )
	{
		return;
	}
	
	const std::vector< Core::Point >& vertices = this->tool_->vertices_state_->get();
	size_t num_of_vertices = vertices.size();
	if ( num_of_vertices < 3 )
	{
		return;
	}
//	double world_x, world_y;
//	int x, y;
//	std::vector< ActionEdgeQuery::VertexCoord > vertices_2d;
//
//	for ( size_t i = 0; i < num_of_vertices; ++i )
//	{
//		volume_slice->project_onto_slice( vertices[ i ], world_x, world_y );
//		volume_slice->world_to_index( world_x, world_y, x, y );
//		vertices_2d.push_back( ActionEdgeQuery::VertexCoord( 
//			static_cast< float >( x ), static_cast< float >( y ), 0 ) );
//	}
//	
//	ActionEdgeQuery::Dispatch( context, this->tool_->target_layer_state_->get(),
//		volume_slice->get_slice_type(), volume_slice->get_slice_number(), save, stop,
//    this->tool_->selectedEdges_state_->get(), vertices_2d );
}

//////////////////////////////////////////////////////////////////////////
// Class EdgeQueryTool
//////////////////////////////////////////////////////////////////////////

EdgeQueryTool::EdgeQueryTool( const std::string& toolid ) :
	SingleTargetTool( Core::VolumeType::MASK_E, toolid ),
	private_( new EdgeQueryToolPrivate )
{
	this->private_->tool_ = this;
	//this->private_->moving_vertex_ = false;

	this->add_state( "vertices", this->vertices_state_ );
  this->add_state( "edges", this->selectedEdges_state_ );
  this->add_state( "save", this->save_state_, false );
  this->add_state( "stop", this->stop_state_, false );
  
	this->add_connection( this->vertices_state_->state_changed_signal_.connect(
    boost::bind( &EdgeQueryToolPrivate::handle_vertices_changed, this->private_ ) ) );
  
	this->add_connection( this->selectedEdges_state_->state_changed_signal_.connect(
    boost::bind( &EdgeQueryToolPrivate::handle_selection_changed, this->private_ ) ) );

//	this->add_connection( this->save_state_->state_changed_signal_.connect(
//    boost::bind( &EdgeQueryToolPrivate::handle_save_changed, this->private_ ) ) );
}

EdgeQueryTool::~EdgeQueryTool()
{
	this->disconnect_all();
}

void EdgeQueryTool::save( Core::ActionContextHandle context )
{
  Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
  this->save_state_->set(true);

  this->private_->execute(context, this->save_state_->get(), this->stop_state_->get());
}

void EdgeQueryTool::stop( Core::ActionContextHandle context )
{
  Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
  this->stop_state_->set(true);
  
  this->private_->execute(context, this->save_state_->get(), this->stop_state_->get());
}

bool EdgeQueryTool::handle_mouse_move( ViewerHandle viewer, 
                        const Core::MouseHistory& mouse_history,
                        int button, int buttons, int modifiers )
{
  this->private_->edgeQuery_.clearHover();

	if ( viewer->is_volume_view() )
	{
		return false;
	}
  
  // May be called from application (slice changed) or interface (mouse move) thread
  //lock_type lock( this->get_mutex() );
  
  Core::VolumeSliceHandle active_slice = viewer->get_active_volume_slice();
  if ( !active_slice )
  {
    return false;
  }
  
  if ( ! this->private_->edgeQuery_.is_valid() )
  {
    return false;
  }

  LineSegment test = this->private_->edgeQuery_.getEdge(0);
  double point_depth = 0;
  if ( ! this->private_->point_in_slice( viewer, test.p1_, point_depth ) )
  {
    return false;
  }  
  
  // Compute the size of a pixel in world space
  double i0, j0, i1, j1;
  viewer->window_to_world( 0, 0, i0, j0 );
  viewer->window_to_world( 1, 1, i1, j1 );
  double pixel_width = Core::Abs( i1 - i0 );
  double epsilon = pixel_width * 4;
  
  // Compute the mouse position in world space
  double mouse_x, mouse_y;
  viewer->window_to_world( mouse_history.current_.x_, mouse_history.current_.y_, mouse_x, mouse_y );
  Core::Point mouse_point( mouse_x, mouse_y, 0 );
  
  Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
  //    const std::vector< Core::Measurement >& measurements = this->tool_->measurements_state_->get();
    
  //bool hovering = false;
  double min_dist = DBL_MAX;

  for( size_t index = 0; index < EDGE_QUERY_SIZE; ++index )
  {
    LineSegment ls = this->private_->edgeQuery_.getEdge(index);
    
    Core::Point p0, p1;
    p0 = ls.p1_;
    p1 = ls.p2_;
    
    // Project points onto slice so that we can detect if we're hovering over a line
    // not in the current slice.
    double p0_i, p0_j;
    active_slice->project_onto_slice( p0, p0_i, p0_j );
    p0 = Core::Point( p0_i, p0_j, 0 );
    double p1_i, p1_j;
    active_slice->project_onto_slice( p1, p1_i, p1_j );
    p1 = Core::Point( p1_i, p1_j, 0 );
    
    // Calculate the shortest distance between the point and the line segment
    double distance = 0.0;
    double l2 = ( p1 - p0 ).length2();  // i.e. |w-v|^2 -  avoid a sqrt
    if( l2 == 0.0 ) // p0 == p1 case
    {
      distance = ( mouse_point - p0 ).length();
    }
    else
    {
      // Consider the line extending the segment, parameterized as v + t (w - v).
      // We find projection of point p onto the line. 
      // It falls where t = [(p-v) . (w-v)] / |w-v|^2
      double t = Dot( mouse_point - p0, p1 - p0 ) / l2;
      if( t < 0.0 || t > 1.0 )
      {
        // Beyond the segment
        continue;
      }
      else
      {
        Core::Point projection = p0 + t * ( p1 - p0 );  // Projection falls on the segment
        distance = ( mouse_point - projection ).length();
      }	
    }
    
    if ( distance > epsilon || distance > min_dist )
    {
      viewer->set_cursor( Core::CursorShape::CROSS_E );
      //this->private_->edgeQuery_.unselect(index);
      this->private_->edgeQuery_.unsetHover(index);
      continue;
    }
    
    // Found the new closest measurement
    //hovering = true;
    min_dist = distance;
    viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );

    //this->private_->edgeQuery_.select(index);
    this->private_->edgeQuery_.hover(index);

    return true;
  }
  
  return false;
}

bool EdgeQueryTool::handle_mouse_press( ViewerHandle viewer, 
									  const Core::MouseHistory& mouse_history, 
									  int button, int buttons, int modifiers )
{
  // ignoring...

//	Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
//	if ( viewer->is_volume_view() )
//	{
//		return false;
//	}
  
  // let default renderer handle mouse press event
	return false;
}

bool EdgeQueryTool::handle_mouse_release( ViewerHandle viewer, 
										const Core::MouseHistory& mouse_history, 
										int button, int buttons, int modifiers )
{
	if ( viewer->is_volume_view() )
	{
		return false;
	}

  // ignore
  if ( button == Core::MouseButton::RIGHT_BUTTON_E || button == Core::MouseButton::MID_BUTTON_E )
  {
    return false;
  }

  // safe, but not efficient...
  this->selectedEdges_state_->clear();

  // TODO: change to bool?
  Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
  for (int i = 0; i < EDGE_QUERY_SIZE; ++i)
  {
    if ( this->private_->edgeQuery_.hovering(i) )
    {
      // toggle selection
      if ( this->private_->edgeQuery_.selected(i) )
      {
        this->private_->edgeQuery_.unselect(i);
        this->selectedEdges_state_->add(0);
      }
      else
      {
        this->private_->edgeQuery_.select(i);
        this->selectedEdges_state_->add(1);
      }
    }
    else
    {
      if ( this->private_->edgeQuery_.selected(i) )
      {
        this->selectedEdges_state_->add(1);
      }
      else
      {
        this->selectedEdges_state_->add(0);
      }
    }
  }

	return false;
}

void EdgeQueryTool::redraw( size_t viewer_id, const Core::Matrix& proj_mat,
	int viewer_width, int viewer_height )
{
	ViewerHandle viewer = ViewerManager::Instance()->get_viewer( viewer_id );
	if ( viewer->is_volume_view() )
	{
		return;
	}

	Core::VolumeSliceType slice_type( Core::VolumeSliceType::AXIAL_E );
	{
		Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );

		if ( viewer->view_mode_state_->get() == Viewer::SAGITTAL_C )
		{
			slice_type = Core::VolumeSliceType::SAGITTAL_E;
		}
		else if ( viewer->view_mode_state_->get() == Viewer::CORONAL_C )
		{
			slice_type = Core::VolumeSliceType::CORONAL_E;
		}
	}
  
  if (! this->private_->edgeQuery_.is_valid() )
  {
    // TODO: debug log
    //std::cerr << "No EdgeQuery available!!!" << std::endl;
    return;
  }
  
  LineSegment ls1 = this->private_->edgeQuery_.getEdge(0);
  LineSegment ls2 = this->private_->edgeQuery_.getEdge(1);
  double point_depth = 0;

  // TODO: see if drawing dashed line in more subdued color is a useful visual cue
  if (! (this->private_->point_in_slice( viewer, ls1.p1_, point_depth ) && this->private_->point_in_slice( viewer, ls1.p2_, point_depth ) ) )
  {
    //std::cerr << "LS 1 not in slice" << std::endl;
    ls1.notInSlice_ = true;
  }
  
  if (! (this->private_->point_in_slice( viewer, ls2.p1_, point_depth ) && this->private_->point_in_slice( viewer, ls2.p2_, point_depth ) ) )
  {
    //std::cerr << "LS 2 not in slice" << std::endl;
    ls2.notInSlice_ = true;
  }
  
  if (ls1.notInSlice_ && ls2.notInSlice_)
  {
    return;
  }
  
  
	glPushAttrib( GL_LINE_BIT | GL_POINT_BIT | GL_TRANSFORM_BIT );
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd( proj_mat.data() );
  
	glPointSize( 6.0f );
	glLineWidth( 4.0f );
	glEnable( GL_LINE_SMOOTH );

	glBegin( GL_POINTS );

  double ls1_p1_iPos, ls1_p1_jPos, ls1_p2_iPos, ls1_p2_jPos, 
         ls2_p1_iPos, ls2_p1_jPos, ls2_p2_iPos, ls2_p2_jPos;

  Core::VolumeSlice::ProjectOntoSlice( slice_type, ls1.p1_, ls1_p1_iPos, ls1_p1_jPos );
  Core::VolumeSlice::ProjectOntoSlice( slice_type, ls1.p2_, ls1_p2_iPos, ls1_p2_jPos );
  Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p1_, ls2_p1_iPos, ls2_p1_jPos );
  Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p2_, ls2_p2_iPos, ls2_p2_jPos );

  if (slice_type == Core::VolumeSliceType::SAGITTAL_E)
  {
    // y-z plane
    ls1.p1_[ 1 ] = ls1_p1_iPos;
    ls1.p1_[ 2 ] = ls1_p1_jPos;
  }
  else if (slice_type == Core::VolumeSliceType::CORONAL_E)
  {
    // x-z plane
    ls1.p1_[ 0 ] = ls1_p1_iPos;
    ls1.p1_[ 2 ] = ls1_p1_jPos;
  }
  else // AXIAL
  {
    // x-y plane
    ls1.p1_[ 0 ] = ls1_p1_iPos;
    ls1.p1_[ 1 ] = ls1_p1_jPos;
  }
  
  if ( ( this->private_->edgeQuery_.selected(0) ) && ( ! this->private_->edgeQuery_.selected(1) ) )
  {
    // draw all points from line 1 and last point of line 2
    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
		glVertex2d( ls1_p1_iPos, ls1_p1_jPos );
		glVertex2d( ls1_p2_iPos, ls1_p2_jPos );

    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls2_p2_iPos, ls2_p2_jPos );
  }

  else if ( ( ! this->private_->edgeQuery_.selected(0) ) && ( this->private_->edgeQuery_.selected(1) ) )
  {
    // draw first point of line 1 and all points from line 2
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls1_p1_iPos, ls1_p1_jPos );

    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
		glVertex2d( ls2_p1_iPos, ls2_p1_jPos );
		glVertex2d( ls2_p2_iPos, ls2_p2_jPos );
  }
  else
  {
    if ( this->private_->edgeQuery_.allSelected() )
    {
      // draw first point of line 1 and all points from line 2
      // white
      glColor3f( 1.0f, 1.0f, 1.0f );
    }
    else if ( this->private_->edgeQuery_.allUnselected() )
    {
      // draw first point of line 1 and all points from line 2
      // cyan
      glColor3f( 0.0f, 1.0f, 1.0f );
    }
    glVertex2d( ls1_p1_iPos, ls1_p1_jPos );
    glVertex2d( ls2_p1_iPos, ls2_p1_jPos );
    glVertex2d( ls2_p2_iPos, ls2_p2_jPos );
  }
  
  glEnd();

	glBegin( GL_LINE_STRIP );

  if ( this->private_->edgeQuery_.selected(0) && ( ! this->private_->edgeQuery_.selected(1) ) )
  {
    // draw all points through all points in line 1 to last point of line 2
    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
		glVertex2d( ls1_p1_iPos, ls1_p1_jPos );
		glVertex2d( ls1_p2_iPos, ls1_p2_jPos );
    
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls2_p2_iPos, ls2_p2_jPos );
  }

  else if ( ( ! this->private_->edgeQuery_.selected(0) ) && this->private_->edgeQuery_.selected(1) )
  {
    // draw line segments from first point of line 1 through all points in line 2
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls1_p1_iPos, ls1_p1_jPos );
    
    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
		glVertex2d( ls2_p1_iPos, ls2_p1_jPos );
		glVertex2d( ls2_p2_iPos, ls2_p2_jPos );
  }
  else
  {
    if ( this->private_->edgeQuery_.allSelected() )
    {
      // draw line segments from first point of line 1 through all points in line 2
      // white
      glColor3f( 1.0f, 1.0f, 1.0f );
    }
    else if ( this->private_->edgeQuery_.allUnselected() )
    {
      // draw line segments from first point of line 1 through all points in line 2
      // cyan
      glColor3f( 0.0f, 1.0f, 1.0f );
    }
    glVertex2d( ls1_p1_iPos, ls1_p1_jPos );
    glVertex2d( ls2_p1_iPos, ls2_p1_jPos );
    glVertex2d( ls2_p2_iPos, ls2_p2_jPos );
  }
  
  glEnd();
  
	glPopMatrix();
	glPopAttrib();
}

bool EdgeQueryTool::has_2d_visual()
{
	return true;
}

} // end namespace Seg3D
