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
	this->add_connection( renderer_->rendering_completed_signal_.connect(
	    boost::bind( &QtRenderWidget::update_texture, this, _1, _2 ) ) );
	this->add_connection( renderer_->redraw_overlay_completed_signal_.connect(
	    boost::bind( &QtRenderWidget::update_overlay_texture, this, _1 ) ) );

	setAutoFillBackground( false );
	setAttribute( Qt::WA_OpaquePaintEvent );
	setAttribute( Qt::WA_NoSystemBackground );
	setMouseTracking( true );
}

QtRenderWidget::~QtRenderWidget()
{
	this->disconnect_all();
}

void QtRenderWidget::update_texture( Utils::TextureHandle texture, bool delay_update )
{
	// if not in the interface thread, post an event to the interface thread
	if ( !Interface::IsInterfaceThread() )
	{
		Interface::PostEvent(
		    boost::bind( &QtRenderWidget::update_texture, this, texture, delay_update ) );
		return;
	}

	SCI_LOG_DEBUG(std::string("QtRenderWidget ") + Utils::ToString(this->viewer_id_)
		+ ": received new texture");
	renderer_texture_ = texture;
	
	if ( !delay_update )
	{
		this->updateGL();
	}
}

void QtRenderWidget::update_overlay_texture( Utils::TextureHandle texture )
{
	// if not in the interface thread, post an event to the interface thread
	if ( !Interface::IsInterfaceThread() )
	{
		Interface::PostEvent(
		    boost::bind( &QtRenderWidget::update_overlay_texture, this, texture ) );
		return;
	}

	SCI_LOG_DEBUG(std::string("QtRenderWidget ") + Utils::ToString(this->viewer_id_)
		+ ": received new overlay texture");
	this->overlay_texture_ = texture;
	
	this->updateGL();
}

void QtRenderWidget::initializeGL()
{
	Utils::RenderResources::Instance()->init_gl();
	glClearColor( 0.5, 0.5, 0.5, 1.0 );
	Utils::Texture::SetActiveTextureUnit( 0 );
	glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

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
	if ( !this->renderer_texture_ || !this->overlay_texture_ )
	{
		return;
	}


	SCI_LOG_DEBUG("Painting texture");

	// draw a window size quad and map the render texture onto it
	QSize view_size = QWidget::size();
	int width = view_size.width();
	int height = view_size.height();

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

	Utils::Texture::lock_type lock( this->renderer_texture_->get_mutex() );
	Utils::Texture::lock_type overlay_tex_lock( this->overlay_texture_->get_mutex() );

	this->renderer_texture_->enable();
	glBegin( GL_QUADS );
	glTexCoord2f( 0.0f, 0.0f );
	glVertex2f( 0.0f, 0.0f );
	glTexCoord2f( 1.0f, 0.0f );
	glVertex2f( width, 0.0f );
	glTexCoord2f( 1.0f, 1.0f );
	glVertex2f( width, height );
	glTexCoord2f( 0.0f, 1.0f );
	glVertex2f( 0.0f, height );
	glEnd();
	this->renderer_texture_->disable();

	// Enable blending to render the overlay texture.
	// NOTE: The overlay texture can NOT be simply rendered using multi-textureing because
	// its color channals already reflect the effect of transparency, and its alpha channal 
	// actually stores the value of "1-alpha"
	glEnable( GL_BLEND );
	glBlendFunc( GL_ONE, GL_SRC_ALPHA );
	this->overlay_texture_->enable();
	glBegin( GL_QUADS );
	glTexCoord2f( 0.0f, 0.0f );
	glVertex2f( 0.0f, 0.0f );
	glTexCoord2f( 1.0f, 0.0f );
	glVertex2f( width, 0.0f );
	glTexCoord2f( 1.0f, 1.0f );
	glVertex2f( width, height );
	glTexCoord2f( 0.0f, 1.0f );
	glVertex2f( 0.0f, height );
	glEnd();
	this->overlay_texture_->disable();
	glDisable( GL_BLEND );

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
		SCI_LOG_DEBUG(std::string("QtRenderWidget ") + Utils::ToString(this->viewer_id_)
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
