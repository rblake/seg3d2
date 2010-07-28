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

// STL includes
#include <iostream>

// Glew includes
#include <GL/glew.h>

// Core includes
#include <Core/Interface/Interface.h>
#include <Core/State/Actions/ActionSet.h>
#include <Core/Utils/Exception.h>
#include <Core/Utils/Log.h>

// QtUtils includes
#include <QtUtils/Utils/QtRenderWidget.h>

namespace QtUtils
{

typedef QPointer< QtRenderWidget > QtRenderWidgetWeakHandle;

class QtRenderWidgetPrivate : public QObject
{
public:
	Core::AbstractViewerHandle viewer_;
	Core::MouseHistory mouse_history_;
	QtRenderWidget* render_widget_;
};


QtRenderWidget::QtRenderWidget( const QGLFormat& format, QWidget* parent, QGLWidget* share, 
		Core::AbstractViewerHandle viewer ) :
	QGLWidget( format, parent, share ),
	private_( new QtRenderWidgetPrivate )
{
	this->private_->viewer_ = viewer;
	this->private_->render_widget_ = this;

	this->setAutoFillBackground( false );
	this->setAttribute( Qt::WA_OpaquePaintEvent );
	this->setAttribute( Qt::WA_NoSystemBackground );
	this->setMouseTracking( true );
	//this->setCursor( Qt::CrossCursor );
	this->setFocusPolicy( Qt::StrongFocus );
}

QtRenderWidget::~QtRenderWidget()
{
	this->disconnect_all();
}

static void UpdateDisplay( QtRenderWidgetWeakHandle qpointer )
{
	if ( !Core::Interface::IsInterfaceThread() )
	{
		Core::Interface::PostEvent( boost::bind( &UpdateDisplay, qpointer ) );
		return;
	}

	if ( !qpointer.isNull() )
	{
		qpointer->updateGL();
	}
}

void QtRenderWidget::initializeGL()
{
	if ( Core::RenderResources::Instance()->valid_render_resources() )
	{
		glClearColor( 0.5, 0.5, 0.5, 1.0 );
		Core::Texture::SetActiveTextureUnit( 0 );
		glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

		QtRenderWidgetWeakHandle qpointer( this );

		this->add_connection( this->private_->viewer_->update_display_signal_.connect(
			boost::bind( &UpdateDisplay, qpointer ) ) );
	}
}

void QtRenderWidget::paintGL()
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	Core::Texture2DHandle texture = this->private_->viewer_->get_texture();
	Core::Texture2DHandle overlay_texture = this->private_->viewer_->get_overlay_texture();

	if ( !texture || !overlay_texture )
	{
		return;
	}

	// draw a window size quad and map the render texture onto it
	QSize view_size = QWidget::size();
	int width = view_size.width();
	int height = view_size.height();

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

	Core::Texture::lock_type texture_lock( texture->get_mutex() );
	Core::Texture::lock_type overlay_texture_lock( overlay_texture->get_mutex() );

	texture->enable();
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
	texture->disable();

	// Enable blending to render the overlay texture.
	// NOTE: The overlay texture can NOT be simply rendered using multi-texturing because
	// its color channels already reflect the effect of transparency, and its alpha channel 
	// actually stores the value of "1-alpha"
	glEnable( GL_BLEND );
	glBlendFunc( GL_ONE, GL_SRC_ALPHA );
	overlay_texture->enable();
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
	overlay_texture->disable();
	glDisable( GL_BLEND );
}

void QtRenderWidget::resizeGL( int width, int height )
{
	glViewport( 0, 0, width, height );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluOrtho2D( 0, width, 0, height );

	this->private_->viewer_->resize( width, height );
}

void QtRenderWidget::mouseMoveEvent( QMouseEvent * event )
{
	event->accept();
	this->private_->mouse_history_.previous_ = this->private_->mouse_history_.current_;
	this->private_->mouse_history_.current_.x_ = event->x();
	this->private_->mouse_history_.current_.y_ = event->y();

	this->private_->viewer_->mouse_move_event( this->private_->mouse_history_, 
		event->button(), event->buttons(), event->modifiers() );
}

void QtRenderWidget::mousePressEvent( QMouseEvent * event )
{
	this->activate_signal_();

	this->private_->mouse_history_.current_.x_ = 
		this->private_->mouse_history_.previous_.x_ = event->x();
	this->private_->mouse_history_.current_.y_ = 
		this->private_->mouse_history_.previous_.y_ = event->y();
	
	if ( event->button() == Qt::LeftButton )
	{
		this->private_->mouse_history_.left_start_.x_ = event->x();
		this->private_->mouse_history_.left_start_.y_ = event->y();
	}
	else if ( event->button() == Qt::RightButton )
	{
		this->private_->mouse_history_.right_start_.x_ = event->x();
		this->private_->mouse_history_.right_start_.y_ = event->y();
	}
	else if ( event->button() == Qt::MidButton )
	{
		this->private_->mouse_history_.mid_start_.x_ = event->x();
		this->private_->mouse_history_.mid_start_.y_ = event->y();
	}

	this->private_->viewer_->mouse_press_event( 
		this->private_->mouse_history_, event->button(), event->buttons(),
	    event->modifiers() );
}

void QtRenderWidget::mouseReleaseEvent( QMouseEvent * event )
{
	this->private_->mouse_history_.previous_ = this->private_->mouse_history_.current_;
	this->private_->mouse_history_.current_.x_ = event->x();
	this->private_->mouse_history_.current_.y_ = event->y();

	this->private_->viewer_->mouse_release_event( this->private_->mouse_history_, event->button(), 
		event->buttons(), event->modifiers() );
}

void QtRenderWidget::wheelEvent( QWheelEvent* event )
{
	this->activate_signal_();

	int delta = Core::RoundUp( event->delta() / 120.0 );
	if ( this->private_->viewer_->wheel_event( delta, event->x(), event->y(), 
		event->buttons(), event->modifiers() ) )
	{
		event->accept();
	}
	else
	{
		event->ignore();
	}
}

void QtRenderWidget::keyPressEvent( QKeyEvent* event )
{
	this->activate_signal_();

	if ( ! ( this->private_->viewer_->key_press_event( event->key(), event->modifiers() ) ) )
	{
		QWidget::keyPressEvent( event );
	}
}

void QtRenderWidget::hideEvent( QHideEvent* event )
{
	if ( !event->spontaneous() )
	{
		Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(), 
			this->private_->viewer_->viewer_visible_state_, false );
	}
}

void QtRenderWidget::showEvent( QShowEvent* event )
{
	if ( !event->spontaneous() )
	{
		Core::ActionSet::Dispatch( Core::Interface::GetWidgetActionContext(), 
			this->private_->viewer_->viewer_visible_state_, true );
	}
}

void QtRenderWidget::enterEvent( QEvent* event )
{
	event->accept();
	setFocus();
	QPoint cursor_pos = this->mapFromGlobal( QCursor::pos() );
	this->private_->viewer_->mouse_enter_event( cursor_pos.x(), cursor_pos.y() );
}

void QtRenderWidget::leaveEvent( QEvent* event )
{
	event->accept();
	clearFocus();
	this->private_->viewer_->mouse_leave_event();
}

} // end namespace QtUtils
