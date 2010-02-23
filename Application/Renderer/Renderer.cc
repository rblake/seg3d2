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

#include <cstdlib>

// application includes
#include <Application/Renderer/Renderer.h>
#include <Application/Renderer/RenderResources.h>
#include <Application/Renderer/UnitCube.h>
#include <Application/ViewerManager/ViewerManager.h>

#include <Utils/EventHandler/DefaultEventHandlerContext.h>
#include <Utils/Geometry/View3D.h>

namespace Seg3D 
{

class RendererEventHandlerContext : public Utils::DefaultEventHandlerContext
{
public:

	RendererEventHandlerContext() : DefaultEventHandlerContext() { }
	virtual ~RendererEventHandlerContext() { }

	virtual void post_event(Utils::EventHandle& event)
	{
		boost::unique_lock<boost::mutex> lock(event_queue_mutex_);

		// discard the previous rendering event
		//if (!event_queue_.empty())
		//{
		//  event_queue_.pop();
		//}

		event_queue_.push(event);
		event_queue_new_event_.notify_all();
	}
};

Renderer::Renderer() : 
	ViewerRenderer(), EventHandler(), 
	active_render_texture_(0), width_(0), 
	height_(0), redraw_needed_(false)
{
}

Renderer::~Renderer() 
{
}

void Renderer::initialize()
{
#if defined(WIN32) || defined(APPLE) || defined(X11_THREADSAFE)
	// NOTE: it is important to postpone the allocation of OpenGL objects to the 
	// rendering thread. If created in a different thread, these objects might not
	// be ready when the rendering thread uses them the first time, which caused
	// the scene to be blank sometimes.
	if (!is_eventhandler_thread())
	{
		if (!RenderResources::Instance()->create_render_context(context_))
		{
			SCI_THROW_EXCEPTION("Failed to create a valid rendering context");
		}
		post_event(boost::bind(&Renderer::initialize, this));
		start_eventhandler();
		return;
	}
#else
	if (!Interface::IsInterfaceThread())
	{
		Interface::PostEvent(boost::bind(&Renderer::initialize, this));
		return;
	}
	if (!RenderResources::Instance()->create_render_context(context_))
	{
		SCI_THROW_EXCEPTION("Failed to create a valid rendering context");
	}
#endif

	SCI_LOG_DEBUG(std::string("Renderer ") + Utils::to_string(this->viewer_id_)
				  + ": initializing");
	
	// Make the GL context current. Since it is the only context in the rendering
	// thread, this call is only needed once.
	this->context_->make_current();

	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);

	// lock the shared render context
	RenderResources::lock_type lock(
		RenderResources::Instance()->shared_context_mutex());

	textures_[0] = TextureHandle(new Texture2D());
	textures_[1] = TextureHandle(new Texture2D());
	depth_buffer_ = RenderBufferHandle(new RenderBuffer());
	frame_buffer_ = FrameBufferObjectHandle(new FrameBufferObject());
	frame_buffer_->attach_render_buffer(depth_buffer_, GL_DEPTH_ATTACHMENT_EXT);
	this->cube_ = UnitCubeHandle(new UnitCube());

	// release the lock
	lock.unlock();

	ViewerManager::Instance()->get_viewer(this->viewer_id_)
		->redraw_signal_.connect(boost::bind(&Renderer::redraw, this));
}

void Renderer::redraw()
{
#if defined(WIN32) || defined(APPLE) || defined(X11_THREADSAFE)
	if (!is_eventhandler_thread())
	{
		boost::unique_lock<boost::recursive_mutex> lock(redraw_needed_mutex_);
		redraw_needed_ = true;
		post_event(boost::bind(&Renderer::redraw, this));
		return;
	}
#else
	if (!Interface::IsInterfaceThread())
	{
		boost::unique_lock<boost::recursive_mutex> lock(redraw_needed_mutex_);
		redraw_needed_ = true;
		Interface::PostEvent(boost::bind(&Renderer::redraw, this));
		return;
	}
#endif

	{
		boost::unique_lock<boost::recursive_mutex> lock(redraw_needed_mutex_);
		redraw_needed_ = false;
	}

	{
		boost::unique_lock<boost::recursive_mutex> lock(redraw_needed_mutex_);
		if (redraw_needed_)
		{
			return;
		}
	}

	// lock the active render texture
	Texture::lock_type texture_lock(textures_[active_render_texture_]->get_mutex());

#if !defined(WIN32) && !defined(APPLE) && !defined(X11_THREADSAFE)
	this->context_->make_current();
#endif

	frame_buffer_->attach_texture(textures_[active_render_texture_]);
	
	SCI_CHECK_OPENGL_ERROR();

	// bind the framebuffer object
	frame_buffer_->enable();
	if ( !frame_buffer_->check_status() )
	{
		return;
	}
	//glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	SCI_CHECK_OPENGL_ERROR();

	glViewport(0, 0, width_, height_); 
	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	SCI_CHECK_OPENGL_ERROR();
	// do some rendering
	// ...
	// ...
	//glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	ViewerHandle viewer = ViewerManager::Instance()->get_viewer(this->viewer_id_);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (viewer->is_volume_view())
	{
		StateEngine::lock_type lock(StateEngine::Instance()->get_mutex());
		Utils::View3D view3d( viewer->volume_view_state_->get() );
		lock.unlock();
		gluPerspective(view3d.fov(), width_ / (1.0 * height_), 0.1, 5.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(view3d.eyep().x(), view3d.eyep().y(), view3d.eyep().z(), 
			view3d.lookat().x(), view3d.lookat().y(), view3d.lookat().z(),
			view3d.up().x(), view3d.up().y(), view3d.up().z());
	}
	else
	{
		StateEngine::lock_type lock(StateEngine::Instance()->get_mutex());
		Utils::View2D view2d( 
			dynamic_cast<StateView2D*>(viewer->get_active_view_state().get())->get() );
		lock.unlock();
		double left, right, top, bottom;
		this->compute_2d_clipping_planes(view2d, left, right, bottom, top);
		gluOrtho2D(left, right, bottom, top);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
	
	SCI_CHECK_OPENGL_ERROR();
	glRotatef(25.0f * (this->viewer_id_+1), 1, 0, 1);
	glScalef(0.5f, 0.5f, 0.5f);
	glTranslatef(-0.5f, -0.5f, -0.5f);
	this->cube_->draw();

	SCI_CHECK_OPENGL_ERROR();
	//glBegin(GL_TRIANGLES);
	//glColor3f(1.0, 0.0, 0.0);
	//glVertex3f(0.5, -0.5, 0);
	//glColor3f(0.0, 1.0, 0.0);
	//glVertex3f(0.0, 0.5, 0);
	//glColor3f(0.0, 0.0, 1.0);
	//glVertex3f(-0.5, -0.5, 0);
	//glEnd();

	//unsigned char* pixels = new unsigned char[(width_)*(height_)*3];
	//glReadPixels(0, 0, width_, height_, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)pixels);
	//unsigned char val = pixels[width_*height_+1];

	glFlush();

	//err = glGetError();

	//if (pixels != NULL)
	//{
	//  delete[] pixels;
	//}

	frame_buffer_->disable(); 

	// release the lock on the active render texture
	texture_lock.unlock();

	glFinish();

	// signal rendering completed
	rendering_completed_signal(textures_[active_render_texture_]);

	// swap render textures 
	active_render_texture_ = (~active_render_texture_)&1;
}

void Renderer::resize(int width, int height)
{
#if defined(WIN32) || defined(APPLE) || defined(X11_THREADSAFE)
	if (!is_eventhandler_thread())
	{
		post_event(boost::bind(&Renderer::resize, this, width, height));
		return;
	}
#else
	if (!Interface::IsInterfaceThread())
	{
		Interface::PostEvent(boost::bind(&Renderer::resize, this, width, height));
		return;
	}
#endif

	if ( width == 0 || height == 0
		|| (width_ == width && height_ == height) )
	{
		return;
	}

	{
		RenderResources::lock_type lock(RenderResources::Instance()->shared_context_mutex());
		textures_[0] = TextureHandle(new Texture2D());
		textures_[1] = TextureHandle(new Texture2D());
		textures_[0]->set_image(width, height, 1, GL_RGBA);
		textures_[1]->set_image(width, height, 1, GL_RGBA);

		depth_buffer_->set_storage(width, height, GL_DEPTH_COMPONENT);
	}

	width_ = width;
	height_ = height;

	redraw();
}

void Renderer::compute_2d_clipping_planes( const Utils::View2D& view2d, 
										  double& left, double& right, double& bottom, double& top )
{
	double dimension = Utils::Min(this->width_, this->height_);
	double clipping_width = this->width_ / dimension / view2d.scalex() * 0.5;
	double clipping_height = this->height_ / dimension / view2d.scaley() * 0.5;
	left = view2d.center().x() - clipping_width;
	right = view2d.center().x() + clipping_width;
	bottom = view2d.center().y() - clipping_height;
	top =view2d.center().y() + clipping_height;
}

} // end namespace Seg3D


