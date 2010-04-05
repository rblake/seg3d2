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

// Glew includes
#include <GL/glew.h>

// Utils includes
#include <Utils/Core/Exception.h>
#include <Utils/Core/Log.h>

// Application includes
#include <Application/Interface/Interface.h>
#include <Application/ViewerManager/ViewerManager.h>

// Interface includes
#include <Interface/QtInterface/QtRenderWidget.h>

namespace Seg3D
{

QtRenderWidget::QtRenderWidget( const QGLFormat& format, QWidget* parent, QtRenderWidget* share ) :
	QGLWidget( format, parent, share )
{
	this->renderer_ = RendererHandle( new Renderer() );
	this->rendering_completed_connection_ = renderer_->rendering_completed_signal_.connect(
	    boost::bind( &QtRenderWidget::update_texture, this, _1 ) );

	setAutoFillBackground( false );
	setAttribute( Qt::WA_OpaquePaintEvent );
	setAttribute( Qt::WA_NoSystemBackground );
	setMouseTracking( true );
}

QtRenderWidget::~QtRenderWidget()
{
	rendering_completed_connection_.disconnect();
}

void QtRenderWidget::update_texture( Utils::TextureHandle texture )
{
	// if not in the interface thread, post an event to the interface thread
	if ( !Interface::IsInterfaceThread() )
	{
		Interface::PostEvent(
		    boost::bind( &QtRenderWidget::update_texture, this, texture ) );
		return;
	}

	SCI_LOG_DEBUG(std::string("QtRenderWidget ") + Utils::to_string(this->viewer_id_)
		+ ": received new texture");
	renderer_texture_ = texture;
	
	updateGL();
}

void QtRenderWidget::initializeGL()
{
	Utils::RenderResources::Instance()->init_gl();
	glClearColor( 0.5, 0.5, 0.5, 1.0 );
	renderer_->initialize();
	// Make sure the GL context of the widget is the current one of this thread,
	// because in the single threaded rendering mode, the renderer will make its own context
	// the current one of the Qt thread.
	this->makeCurrent();
}

void QtRenderWidget::paintGL()
{
	SCI_LOG_DEBUG("Start of paintGL");
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	if ( !( this->renderer_texture_ ) )
	{
		return;
	}

	Utils::Texture::lock_type lock( renderer_texture_->get_mutex() );

	SCI_LOG_DEBUG("Painting texture");

	// draw a window size quad and map the render texture onto it
	QSize view_size = QWidget::size();
	int width = view_size.width();
	int height = view_size.height();

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	renderer_texture_->enable();
	glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
	//GLenum err = glGetError();
	//const GLubyte* err_str = gluErrorString(err);
	glBegin( GL_QUADS );
	glColor3f( 0.5f, 0.5f, 1.f );
	glTexCoord2f( 0.0f, 0.0f );
	glVertex2f( 0.0f, 0.0f );
	glTexCoord2f( 1.0f, 0.0f );
	glVertex2f( width, 0.0f );
	glTexCoord2f( 1.0f, 1.0f );
	glVertex2f( width, height );
	glTexCoord2f( 0.0f, 1.0f );
	glVertex2f( 0.0f, height );
	glEnd();
	renderer_texture_->disable();
}

void QtRenderWidget::resizeGL( int width, int height )
{
	glViewport( 0, 0, width, height );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluOrtho2D( 0, width, 0, height );

	this->viewer_->resize( width, height );
	if ( renderer_.get() )
	{
		SCI_LOG_DEBUG(std::string("QtRenderWidget ") + Utils::to_string(this->viewer_id_)
			+ ": sending resize event to renderer");
		renderer_->resize( width, height );
		// Make sure the GL context of the widget is the current one of this thread,
		// because in the single threaded rendering mode, the renderer will make its own context
		// the current one of the Qt thread.
		this->makeCurrent();
	}
}

void QtRenderWidget::mouseMoveEvent( QMouseEvent * event )
{
	mouse_history_.previous = mouse_history_.current;
	mouse_history_.current.x = event->x();
	mouse_history_.current.y = event->y();

	viewer_->mouse_move_event( this->mouse_history_, event->button(), event->buttons(),
	    event->modifiers() );
}

void QtRenderWidget::mousePressEvent( QMouseEvent * event )
{
	mouse_history_.current.x = mouse_history_.previous.x = event->x();
	mouse_history_.current.y = mouse_history_.previous.y = event->y();
	if ( event->button() == Qt::LeftButton )
	{
		mouse_history_.left_start.x = event->x();
		mouse_history_.left_start.y = event->y();
	}
	else if ( event->button() == Qt::RightButton )
	{
		mouse_history_.right_start.x = event->x();
		mouse_history_.right_start.y = event->y();
	}
	else if ( event->button() == Qt::MidButton )
	{
		mouse_history_.mid_start.x = event->x();
		mouse_history_.mid_start.y = event->y();
	}

	viewer_->mouse_press_event( this->mouse_history_, event->button(), event->buttons(),
	    event->modifiers() );
}

void QtRenderWidget::mouseReleaseEvent( QMouseEvent * event )
{
	mouse_history_.previous = mouse_history_.current;
	mouse_history_.current.x = event->x();
	mouse_history_.current.y = event->y();

	viewer_->mouse_release_event( this->mouse_history_, event->button(), event->buttons(),
	    event->modifiers() );
}

void QtRenderWidget::wheelEvent( QWheelEvent* event )
{
	int delta = Utils::RoundUp( event->delta() / 120.0 );
	if ( this->viewer_->wheel_event( delta, event->x(), event->y(), 
		event->buttons(), event->modifiers() ) )
	{
		event->accept();
	}
	else
	{
		event->ignore();
	}
}

void QtRenderWidget::hideEvent( QHideEvent* event )
{
	if ( !event->spontaneous() )
	{
		this->renderer_->deactivate();
	}
}

void QtRenderWidget::showEvent( QShowEvent* event )
{
	if ( !event->spontaneous() )
	{
		this->renderer_->activate();
	}
}

void QtRenderWidget::set_viewer_id( size_t viewer_id )
{
	this->viewer_id_ = viewer_id;
	viewer_ = ViewerManager::Instance()->get_viewer( viewer_id );
	renderer_->set_viewer_id( viewer_id );
}

} // end namespace Seg3D
