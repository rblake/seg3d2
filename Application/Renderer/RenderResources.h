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

#ifndef APPLICATION_RENDERER_RENDERRESOURCES_H
#define APPLICATION_RENDERER_RENDERRESOURCES_H

// Boost includes
#include <boost/thread/mutex.hpp>
#include <boost/thread/recursive_mutex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

// Utils includes
#include <Utils/Core/Exception.h>
#include <Utils/Core/Log.h>
#include <Utils/Singleton/Singleton.h>

// Application includes
#include <Application/Renderer/RenderContext.h>
#include <Application/Renderer/RenderResourcesContext.h>

namespace Seg3D {

class RenderResources : public boost::noncopyable {

// -- constructor --
  public:
    RenderResources();
    
// -- context handling --    
  public:  
    // CREATE_RENDER_CONTEXT:
    // Generate a render context for one of the viewers
    bool create_render_context(RenderContextHandle& context);
    
    // INSTALL_RESOURCES_CONTEXT:
    // Install a context from the UI system in the resource, so this
    // class can generate virtual OpenGL contexts. 
    void install_resources_context(RenderResourcesContextHandle resources_context);
    
    // VALID_RENDER_RESOURCES:
    // Check whether valid render resources have been installed
    bool valid_render_resources();
    
  private:
    
    // A Handle to resource that generated the contexts
    RenderResourcesContextHandle resources_context_;
    
// -- render locks --
  public:

    // LOCK_SHARED_CONTEXT:
    void lock_shared_context() { shared_context_mutex_.lock(); }

    // LOCK_SHARED_CONTEXT:
    void unlock_shared_context() { shared_context_mutex_.unlock(); }

    // SHARED_CONTEXT_MUTEX:
    // reference to the mutex that protects the common pool of textures
    boost::recursive_mutex& shared_context_mutex() { return shared_context_mutex_; }

  private:
    // We need a lock that protects against multiple threads reserving
    // OpenGL resources at the same time
    boost::recursive_mutex shared_context_mutex_;

// -- Singleton interface --
  public:
    
    static RenderResources* Instance() { return instance_.instance(); }
    
  private:
    static Utils::Singleton<RenderResources> instance_;
};

} // end namespace Seg3D

#endif
