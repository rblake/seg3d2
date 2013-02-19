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
  LineSegment(const Core::Point& p1, const Core::Point& p2)
    : p1_(p1), p2_(p2), selected_(false), hovering_(false)
  {}  
  
  Core::Point p1_;
  Core::Point p2_;
  
  bool selected_;
  bool hovering_;
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

  
  bool unselected() const
  {
    return (! lines_[0].selected_) && (! lines_[1].selected_);
  }
  
  void create_edgeQuery(const Core::Point& p1, const Core::Point& p2, const Core::Point& p3)
  {
    lines_.reserve(EDGE_QUERY_SIZE);

    LineSegment e1(p1, p2);
    LineSegment e2(p2, p3);
    lines_.push_back(e1);
    lines_.push_back(e2);
  }
  
  int selectedEdge() const
  {
    if (lines_.size() < EDGE_QUERY_SIZE)
      return -1;

    for (size_t i = 0; i < EDGE_QUERY_SIZE; ++i)
    {
      if ( lines_[i].selected_ ) return i;
    }

    return -1;
  }

  void unselectEdge(int index)
  {
    if (index > EDGE_QUERY_SIZE - 1)
    {
      return;
    }
    lines_[index].selected_ = false;
  }
  
  void selectEdge(int index)
  {
    if (index > EDGE_QUERY_SIZE - 1)
    {
      return;
    }

    for (int i = 0; i < EDGE_QUERY_SIZE; ++i)
    {
      if (i == index)
        lines_[i].selected_ = true;
      else
        lines_[i].selected_ = false;
    }
  }
  
  LineSegment getEdge(int index) const
  {
    // TODO: should have error checking here...exception?
    // this implementation assumes validity already checked

    return lines_[index];
  }
  
private:
  std::vector<LineSegment> lines_;
};


//////////////////////////////////////////////////////////////////////////
// Class EdgeQueryToolPrivate
//////////////////////////////////////////////////////////////////////////

class EdgeQueryToolPrivate
{
public:
	void handle_vertices_changed();
  void handle_selected_edge_changed();
//  void handle_save_changed();

	void execute( Core::ActionContextHandle context, bool erase, 
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
  
  if ( point_depth == slice_depth ) 
  {
    return true;
  }
  return false;
}
  
void EdgeQueryToolPrivate::handle_vertices_changed()
{
  std::cout << "EdgeQueryToolPrivate::handle_vertices_changed()" << std::endl;

  //Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );

  std::vector< Core::Point > points = this->tool_->vertices_state_->get();
  if (points.size() >= 3)
  {
    this->edgeQuery_.create_edgeQuery(points[0], points[1], points[2]);
 
    Core::ActionSet::Dispatch( Core::Interface::GetMouseActionContext(),
                              this->tool_->vertices_state_, points );
   
    
    ViewerManager::Instance()->update_2d_viewers_overlay();
    
//    {
//      Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
//      this->save_state_->set(false);
//    }

  }
}

void EdgeQueryToolPrivate::handle_selected_edge_changed()
{
  std::cout << "EdgeQueryToolPrivate::handle_selected_edge_changed(): " << this->tool_->selectedEdge_state_->get() << std::endl;
	ViewerManager::Instance()->update_2d_viewers_overlay();
}

//void EdgeQueryToolPrivate::handle_save_changed()
//{
//  std::cout << "EdgeQueryToolPrivate::handle_save_changed()" << std::endl;
//}

//bool EdgeQueryToolPrivate::find_vertex( ViewerHandle viewer, int x, int y, int& index )
//{
//	// Step 1. Compute the size of a pixel in world space
//	double x0, y0, x1, y1;
//	viewer->window_to_world( 0, 0, x0, y0 );
//	viewer->window_to_world( 1, 1, x1, y1 );
//	double pixel_width = Core::Abs( x1 - x0 );
//	double pixel_height = Core::Abs( y1 - y0 );
//
//	// Step 2. Compute the mouse position in world space
//	double world_x, world_y;
//	viewer->window_to_world( x, y, world_x, world_y );
//
//	// Step 3. Search for the first vertex that's within 2 pixels of current mouse position
//	double range_x = pixel_width * 4;
//	double range_y = pixel_height * 4;
//	std::vector< Core::Point > vertices = this->tool_->vertices_state_->get();
//	Core::VolumeSliceType slice_type( Core::VolumeSliceType::AXIAL_E );
//	if ( viewer->view_mode_state_->get() == Viewer::CORONAL_C )
//	{
//		slice_type = Core::VolumeSliceType::CORONAL_E;
//	}
//	else if ( viewer->view_mode_state_->get() == Viewer::SAGITTAL_C )
//	{
//		slice_type = Core::VolumeSliceType::SAGITTAL_E;
//	}
//	
//	for ( size_t i = 0; i < vertices.size(); ++i )
//	{
//		double pt_x, pt_y;
//		Core::VolumeSlice::ProjectOntoSlice( slice_type, vertices[ i ], pt_x, pt_y );
//		if ( Core::Abs( pt_x - world_x ) <= range_x &&
//			Core::Abs( pt_y - world_y ) <= range_y )
//		{
//			index = static_cast< int >( i );
//			return true;
//		}
//	}
//
//	index = -1;
//	return false;
//}
//
//bool EdgeQueryToolPrivate::find_closest_vertex( ViewerHandle viewer, int x, int y, int& index )
//{
//	// Step 1. Compute the mouse position in world space
//	double world_x, world_y;
//	viewer->window_to_world( x, y, world_x, world_y );
//	
//	// Step 2. Search for the closest vertex to the current mouse position
//	std::vector< Core::Point > vertices = this->tool_->vertices_state_->get();
//	Core::VolumeSliceType slice_type( Core::VolumeSliceType::AXIAL_E );
//	if ( viewer->view_mode_state_->get() == Viewer::CORONAL_C )
//	{
//		slice_type = Core::VolumeSliceType::CORONAL_E;
//	}
//	else if ( viewer->view_mode_state_->get() == Viewer::SAGITTAL_C )
//	{
//		slice_type = Core::VolumeSliceType::SAGITTAL_E;
//	}
//	
//	int closest_index = -1;
//	double min_dist;
//	for ( size_t i = 0; i < vertices.size(); ++i )
//	{
//		double pt_x, pt_y;
//		Core::VolumeSlice::ProjectOntoSlice( slice_type, vertices[ i ], pt_x, pt_y );
//		double dist_x = pt_x - world_x;
//		double dist_y = pt_y - world_y;
//		double distance = dist_x * dist_x + dist_y * dist_y;
//		if ( i == 0 )
//		{
//			closest_index = 0;
//			min_dist = distance;
//		}
//		else if ( distance < min_dist )
//		{
//			closest_index = static_cast< int >( i );
//			min_dist = distance;
//		}
//	}
//	
//	index = closest_index;
//	return index >= 0;
//}
	
void EdgeQueryToolPrivate::execute( Core::ActionContextHandle context, 
								  bool save, ViewerHandle viewer )
{
  std::cerr << "EdgeQueryToolPrivate::execute" << std::endl;

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
	double world_x, world_y;
	int x, y;
	std::vector< ActionEdgeQuery::VertexCoord > vertices_2d;

  // need to actually set
  int edge = -1;

	for ( size_t i = 0; i < num_of_vertices; ++i )
	{
		volume_slice->project_onto_slice( vertices[ i ], world_x, world_y );
		volume_slice->world_to_index( world_x, world_y, x, y );
		vertices_2d.push_back( ActionEdgeQuery::VertexCoord( 
			static_cast< float >( x ), static_cast< float >( y ), 0 ) );
	}
	
	ActionEdgeQuery::Dispatch( context, this->tool_->target_layer_state_->get(),
		volume_slice->get_slice_type(), volume_slice->get_slice_number(), save, edge, vertices_2d );
  
  // clear selection?
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
  this->add_state( "edge", this->selectedEdge_state_, -1 );
  this->add_state( "save", this->save_state_, false );
  
	this->add_connection( this->vertices_state_->state_changed_signal_.connect(
    boost::bind( &EdgeQueryToolPrivate::handle_vertices_changed, this->private_ ) ) );
  
	this->add_connection( this->selectedEdge_state_->state_changed_signal_.connect(
    boost::bind( &EdgeQueryToolPrivate::handle_selected_edge_changed, this->private_ ) ) );

//	this->add_connection( this->selectedEdge_state_->state_changed_signal_.connect(
//    boost::bind( &EdgeQueryToolPrivate::handle_save_changed, this->private_ ) ) );
}

EdgeQueryTool::~EdgeQueryTool()
{
	this->disconnect_all();
}

void EdgeQueryTool::save( Core::ActionContextHandle context )
{
  Core::StateEngine::lock_type lock( Core::StateEngine::GetMutex() );
std::cerr << "save=" << this->save_state_->export_to_string() << std::endl;
  this->save_state_->set(true);

  this->private_->execute(context, true);
}


bool EdgeQueryTool::handle_mouse_move( ViewerHandle viewer, 
                        const Core::MouseHistory& mouse_history,
                        int button, int buttons, int modifiers )
{
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
  double x0, y0, x1, y1;
  viewer->window_to_world( 0, 0, x0, y0 );
  viewer->window_to_world( 1, 1, x1, y1 );
  double pixel_width = Core::Abs( x1 - x0 );
  double epsilon = pixel_width * 4;
  
  // Compute the mouse position in world space
  double mouse_x, mouse_y;
  viewer->window_to_world( mouse_history.current_.x_, mouse_history.current_.y_, mouse_x, mouse_y );
  Core::Point mouse_point( mouse_x, mouse_y, 0 );
  
  Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
  //    const std::vector< Core::Measurement >& measurements = this->tool_->measurements_state_->get();
    
  bool hovering = false;
  double min_dist = DBL_MAX;

  for( size_t i = 0; i < EDGE_QUERY_SIZE; ++i )
  {
    LineSegment ls = this->private_->edgeQuery_.getEdge(i);
    
    Core::Point p0, p1;
    p0 = ls.p1_;
    p1 = ls.p2_;
    
    // Project points onto slice so that we can detect if we're hovering over a line
    // not in the current slice.
    double p0_x, p0_y;
    active_slice->project_onto_slice( p0, p0_x, p0_y );
    p0 = Core::Point( p0_x, p0_y, 0 );
    double p1_x, p1_y;
    active_slice->project_onto_slice( p1, p1_x, p1_y );
    p1 = Core::Point( p1_x, p1_y, 0 );
    
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
      this->private_->edgeQuery_.unselectEdge(i);
//std::cerr << "EdgeQueryTool::handle_mouse_move: selectedEdge=" << this->private_->edgeQuery_.selectedEdge() << std::endl;
      continue;
    }
    
    // Found the new closest measurement
    hovering = true;
    min_dist = distance;
    viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );
    this->private_->edgeQuery_.selectEdge(i);
//std::cerr << "EdgeQueryTool::handle_mouse_move: hovering over " << i << ", selectedEdge=" << this->private_->edgeQuery_.selectedEdge() << std::endl;
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
	if ( viewer->is_volume_view() )
	{
		return false;
	}
//  
//  if (button == Core::MouseButton::LEFT_BUTTON_E)
//  {
//  	// May be called from application (slice changed) or interface (mouse move) thread
//    //lock_type lock( this->get_mutex() );
//    
//    Core::VolumeSliceHandle active_slice = viewer->get_active_volume_slice();
//    if( !active_slice )
//    {
//      return false;
//    }
//    
//    // Compute the size of a pixel in world space
//    double x0, y0, x1, y1;
//    viewer->window_to_world( 0, 0, x0, y0 );
//    viewer->window_to_world( 1, 1, x1, y1 );
//    double pixel_width = Core::Abs( x1 - x0 );
//    double epsilon = pixel_width * 4;
//    
//    // Compute the mouse position in world space
//    double mouse_x, mouse_y;
//    viewer->window_to_world( mouse_history.current_.x_, mouse_history.current_.y_, mouse_x, mouse_y );
//    Core::Point mouse_point( mouse_x, mouse_y, 0 );
//    
//    Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
////    const std::vector< Core::Measurement >& measurements = this->tool_->measurements_state_->get();
//    
//    std::vector< Core::Point > points = this->vertices_state_->get();
//    
//    
//    bool hovering = false;
//    double min_dist = DBL_MAX;
////    for( size_t m_idx = 0; m_idx < measurements.size(); m_idx++ )
//    for( size_t i = 0; i < 1; ++i )
//    {
////      const Core::Measurement& m = measurements[ m_idx ];
//      
//      // Get both measurement points
//      Core::Point p0, p1;
////      m.get_point( 0, p0 );
////      m.get_point( 1, p1 );
//
//      
//      p0 = points[0];
//      p1 = points[1];
//      
//      // Project points onto slice so that we can detect if we're hovering over a measurement
//      // not in the current slice.
//      double p0_x, p0_y;
//      active_slice->project_onto_slice( p0, p0_x, p0_y );
//      p0 = Core::Point( p0_x, p0_y, 0 );
//      double p1_x, p1_y;
//      active_slice->project_onto_slice( p1, p1_x, p1_y );
//      p1 = Core::Point( p1_x, p1_y, 0 );
//      
//      // Calculate the shortest distance between the point and the line segment
//      double distance = 0.0;
//      double l2 = ( p1 - p0 ).length2();  // i.e. |w-v|^2 -  avoid a sqrt
//      if( l2 == 0.0 ) // p0 == p1 case
//      {
//        distance = ( mouse_point - p0 ).length();
//      }
//      else
//      {
//        // Consider the line extending the segment, parameterized as v + t (w - v).
//        // We find projection of point p onto the line. 
//        // It falls where t = [(p-v) . (w-v)] / |w-v|^2
//        double t = Dot( mouse_point - p0, p1 - p0 ) / l2;
//        if( t < 0.0 || t > 1.0 )
//        {
//          // Beyond the segment
//          continue;
//        }
//        else
//        {
//          Core::Point projection = p0 + t * ( p1 - p0 );  // Projection falls on the segment
//          distance = ( mouse_point - projection ).length();
//        }	
//      }
//      
//      if( distance > epsilon || distance > min_dist )
//      {
//        continue;
//      }
//      
//      // Found the new closest measurement
//      hovering = true;
//      min_dist = distance;
//      //this->hover_measurement_.index_ = static_cast< int >( m_idx );
//      viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );
//      //this->hover_measurement_.hover_object_ = MeasurementHoverObject::LINE_E;
//    }
//    
//    return hovering;
//  }
  
  
  //  if ( button == Core::MouseButton::LEFT_BUTTON_E && this->private_->vertex_index_ != -1 )
//  {
//    return true;
//  }
//
//
//	if ( button == Core::MouseButton::LEFT_BUTTON_E &&
//		( modifiers == Core::KeyModifier::NO_MODIFIER_E ||
//		modifiers == Core::KeyModifier::SHIFT_MODIFIER_E ) &&
//		this->private_->vertex_index_ != -1 )
//	{
//		this->private_->moving_vertex_ = true;
//		viewer->set_cursor( Core::CursorShape::CLOSED_HAND_E );
//		return true;
//	}
//	else if ( button == Core::MouseButton::MID_BUTTON_E &&
//			 ( modifiers == Core::KeyModifier::NO_MODIFIER_E ||
//			  modifiers == Core::KeyModifier::SHIFT_MODIFIER_E ) )
//	{
//		if ( this->private_->find_closest_vertex( viewer, mouse_history.current_.x_, 
//			mouse_history.current_.y_, this->private_->vertex_index_ ) )
//		{
//			this->private_->moving_vertex_ = true;
//			viewer->set_cursor( Core::CursorShape::CLOSED_HAND_E );
//			return true;			
//		}
//	}
//	else if ( !( modifiers & Core::KeyModifier::SHIFT_MODIFIER_E ) &&
//		button == Core::MouseButton::LEFT_BUTTON_E )
//	{
//		Core::VolumeSliceHandle active_slice = viewer->get_active_volume_slice();
//		if ( active_slice && !active_slice->out_of_boundary() )
//		{
//			double world_x, world_y;
//			viewer->window_to_world( mouse_history.current_.x_, 
//				mouse_history.current_.y_, world_x, world_y );
//			Core::Point pt;
//			active_slice->get_world_coord( world_x, world_y,  pt );
//			
//			double dmin = DBL_MAX;
//			double proj_min = DBL_MAX;
//			std::vector<Core::Point> points = this->vertices_state_->get();
//			
//			size_t idx = 0;
//			for ( size_t j = 0; j < points.size(); j++ )
//			{
//				size_t k = j + 1;
//				if ( k ==  points.size() ) k = 0;
//				
//				Core::Vector edge_dir = points[ j ] - points[ k ];
//				double edge_length = edge_dir.normalize();
//				double alpha = Dot( points[ j ] - pt, points[ j ] - points[ k ] )/
//					( edge_length * edge_length );
//					
//				double dist = 0.0;
//				double proj_len = 0.0;
//				if ( alpha < 0.0 ) 
//				{
//					Core::Vector dir = points[ j ] - pt;
//					dist = dir.length2();
//					proj_len = Core::Abs( Dot( edge_dir, dir ) );
//				}
//				else if ( alpha > 1.0 )
//				{
//					Core::Vector dir = points[ k ] - pt;
//					dist = dir.length2();
//					proj_len = Core::Abs( Dot( edge_dir, dir ) );
//				}
//				else 
//				{
//					dist = ( ( points[ j ] - pt ) - alpha * ( points[ j ] - points[ k ] ) ).length2();
//				}
//				
//				if ( dist < dmin ||
//					( dist == dmin && proj_len < proj_min ) )
//				{
//					dmin = dist;
//					proj_min = proj_len;
//					idx = k;
//				}
//			}
//			points.insert( points.begin() + idx, pt );
//			
//			Core::ActionSet::Dispatch( Core::Interface::GetMouseActionContext(),
//				this->vertices_state_, points );
//
//			// Set to "hovered over" state since the mouse is hovering over the new point
//			viewer->set_cursor( Core::CursorShape::OPEN_HAND_E );
//			this->private_->vertex_index_ = static_cast< int >( idx );
//
//			return true;
//		}
//	}
//	else if ( modifiers == Core::KeyModifier::NO_MODIFIER_E &&
//		button == Core::MouseButton::LEFT_BUTTON_E )
//	{
//		Core::VolumeSliceHandle active_slice = viewer->get_active_volume_slice();
//		if ( active_slice && !active_slice->out_of_boundary() )
//		{
//			double world_x, world_y;
//			viewer->window_to_world( mouse_history.current_.x_, 
//				mouse_history.current_.y_, world_x, world_y );
//			Core::Point pt;
//			active_slice->get_world_coord( world_x, world_y,  pt );
//			Core::ActionAdd::Dispatch( Core::Interface::GetMouseActionContext(),
//				this->vertices_state_, pt );
//			return true;
//		}
//	}
//	else if ( modifiers == Core::KeyModifier::NO_MODIFIER_E &&
//		button == Core::MouseButton::RIGHT_BUTTON_E )
//	{
//		if ( this->private_->vertex_index_ != -1 )
//		{
//			Core::Point pt = this->vertices_state_->get()[ this->private_->vertex_index_ ];
//			Core::ActionRemove::Dispatch( Core::Interface::GetMouseActionContext(),
//				this->vertices_state_, pt );
//
//			// Set to "not hovered over" state since the point no longer exists
//			viewer->set_cursor( Core::CursorShape::CROSS_E );
//			this->private_->vertex_index_ = -1;
//
//			return true;
//		}		
//	}

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
  
  //std::cerr << "EdgeQueryTool::handle_mouse_release: selected edge=" << this->private_->edgeQuery_.selectedEdge() << std::endl;

  Core::StateEngine::lock_type state_lock( Core::StateEngine::GetMutex() );
  this->selectedEdge_state_->set(this->private_->edgeQuery_.selectedEdge());

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
    //std::cerr << "No EdgeQuery available!!!" << std::endl;
    return;
  }

  // assumes all points are in the same slice
  // TODO: will this always be true?
  LineSegment ls = this->private_->edgeQuery_.getEdge(0);
  double point_depth = 0;
  if (! this->private_->point_in_slice( viewer, ls.p1_, point_depth ) )
  {
    return;
  }
  
//  Core::VolumeSliceHandle volume_slice = viewer->get_active_volume_slice();
//  if ( ! volume_slice )
//  {
//    return;
//  }
//std::cerr << "volume slice number=" << volume_slice->get_slice_number() << std::endl;
//std::cerr << "point z coord=" << ls.p1_.z() << std::endl;
  
	glPushAttrib( GL_LINE_BIT | GL_POINT_BIT | GL_TRANSFORM_BIT );
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd( proj_mat.data() );

	glPointSize( 6.0f );
	glLineWidth( 4.0f );
  // cyan
	glColor3f( 0.0f, 1.0f, 1.0f );
	glEnable( GL_LINE_SMOOTH );

  LineSegment ls1 = this->private_->edgeQuery_.getEdge(0);
  LineSegment ls2 = this->private_->edgeQuery_.getEdge(1);

	glBegin( GL_POINTS );

  double x_pos, y_pos;

  if ( this->private_->edgeQuery_.selectedEdge() == 0 )
  {
    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls1.p1_, x_pos, y_pos );
		ls1.p1_[ 0 ] = x_pos;
		ls1.p1_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );

		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls1.p2_, x_pos, y_pos );
		ls1.p2_[ 0 ] = x_pos;
		ls1.p2_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );

    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p2_, x_pos, y_pos );
		ls2.p2_[ 0 ] = x_pos;
		ls2.p2_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );
  }
  else if ( this->private_->edgeQuery_.selectedEdge() == 1 )
  {
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls1.p1_, x_pos, y_pos );
		ls1.p1_[ 0 ] = x_pos;
		ls1.p1_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );

    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p1_, x_pos, y_pos );
		ls2.p1_[ 0 ] = x_pos;
		ls2.p1_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );

		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p2_, x_pos, y_pos );
		ls2.p2_[ 0 ] = x_pos;
		ls2.p2_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );
  }
  else
  {
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls1.p1_, x_pos, y_pos );
		ls1.p1_[ 0 ] = x_pos;
		ls1.p1_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );

		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p1_, x_pos, y_pos );
		ls2.p1_[ 0 ] = x_pos;
		ls2.p1_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );
    
		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls2.p2_, x_pos, y_pos );
		ls2.p2_[ 0 ] = x_pos;
		ls2.p2_[ 1 ] = y_pos;
		glVertex2d( x_pos, y_pos );
  }

//  for ( size_t i = 0; i < EDGE_QUERY_SIZE; ++i)
//  {
//    LineSegment ls = this->private_->edgeQuery_.getEdge(i);
//		double x_pos, y_pos;
//
//		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls.p1_, x_pos, y_pos );
//		ls.p1_[ 0 ] = x_pos;
//		ls.p1_[ 1 ] = y_pos;
//		glVertex2d( x_pos, y_pos );
//
//		Core::VolumeSlice::ProjectOntoSlice( slice_type, ls.p2_, x_pos, y_pos );
//		ls.p2_[ 0 ] = x_pos;
//		ls.p2_[ 1 ] = y_pos;
//		glVertex2d( x_pos, y_pos );
//  }
  
  glEnd();

	glBegin( GL_LINE_STRIP );

  //LineSegment ls1 = this->private_->edgeQuery_.getEdge(0);
  //LineSegment ls2 = this->private_->edgeQuery_.getEdge(1);

  if ( this->private_->edgeQuery_.selectedEdge() == 0 )
  {
    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
    glVertex2d( ls1.p1_.x(), ls1.p1_.y() );
    glVertex2d( ls1.p2_.x(), ls1.p2_.y() );
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls2.p2_.x(), ls2.p2_.y() );    
  }
  else if ( this->private_->edgeQuery_.selectedEdge() == 1 )
  {
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls1.p1_.x(), ls1.p1_.y() );
    // white
    glColor3f( 1.0f, 1.0f, 1.0f );
    glVertex2d( ls2.p1_.x(), ls2.p1_.y() );    
    glVertex2d( ls2.p2_.x(), ls2.p2_.y() );    
  }
  else
  {
    // cyan
    glColor3f( 0.0f, 1.0f, 1.0f );
    glVertex2d( ls1.p1_.x(), ls1.p1_.y() );
    glVertex2d( ls1.p2_.x(), ls1.p2_.y() );
    glVertex2d( ls2.p2_.x(), ls2.p2_.y() );
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
