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
// Core includes
#include <Core/Action/ActionDispatcher.h>

// Application includes
#include <Application/ProjectManager/AutoSave.h>
#include <Application/ProjectManager/ProjectManager.h>
#include <Application/ProjectManager/Actions/ActionSaveSession.h>
#include <Application/PreferencesManager/PreferencesManager.h>


namespace Seg3D
{

CORE_SINGLETON_IMPLEMENTATION( AutoSave );

AutoSave::AutoSave()
{
}

AutoSave::~AutoSave()
{
}

void AutoSave::start()
{
	add_connection( PreferencesManager::Instance()->auto_save_time_state_->state_changed_signal_.connect(
		boost::bind( &AutoSave::recompute_auto_save, this ) ) );
	add_connection( PreferencesManager::Instance()->auto_save_state_->state_changed_signal_.connect(
		boost::bind( &AutoSave::recompute_auto_save, this ) ) );
	add_connection( PreferencesManager::Instance()->smart_save_state_->state_changed_signal_.connect(
		boost::bind( &AutoSave::recompute_auto_save, this ) ) );

	this->auto_save_thread_ = boost::thread( boost::bind( &AutoSave::run, this ) );
}

void AutoSave::run()
{
	lock_type lock( this->get_mutex() );
	while( true )
	{	
		double timeout;
		if(	this->needs_auto_save() )
		{
			if( PreferencesManager::Instance()->smart_save_state_->get() == true )
			{
				while( true )
				{
					if( !Core::ActionDispatcher::Instance()->is_busy() )
					{
						boost::posix_time::ptime last_action_completed = 
							Core::ActionDispatcher::Instance()->last_action_completed();
						boost::posix_time::ptime current_time = 
							boost::posix_time::second_clock::local_time();

						boost::posix_time::time_duration duration = current_time - last_action_completed;
						double time_since_last_action = static_cast< double >( duration.total_milliseconds() ) * 0.001;

						if( time_since_last_action > 5  )
						{
							break;
						}
					}
					boost::posix_time::milliseconds action_wait_time( static_cast< int >( 5000.0 ) );
					recompute_auto_save_.timed_wait( lock, action_wait_time );
				}
			}
			this->do_auto_save();
			Core::StateEngine::lock_type state_engine_lock( Core::StateEngine::GetMutex() );
			timeout = PreferencesManager::Instance()->auto_save_time_state_->get() * 60;
		}
		else
		{
			timeout = this->compute_timeout();
		}

		boost::posix_time::milliseconds wait_time( static_cast< int >( 1000.0 * timeout ) );
		recompute_auto_save_.timed_wait( lock, wait_time );	
	}
}

void AutoSave::recompute_auto_save()
{
	lock_type lock( this->get_mutex() );
	this->recompute_auto_save_.notify_all();
}

double AutoSave::compute_timeout()
{
	Core::StateEngine::lock_type state_engine_lock( Core::StateEngine::GetMutex() );
	double time_remaining = PreferencesManager::Instance()->auto_save_time_state_->get() * 60;

	time_remaining = time_remaining - ProjectManager::Instance()->get_time_since_last_saved_session();

	return time_remaining;
}

bool AutoSave::needs_auto_save()
{
	Core::StateEngine::lock_type state_engine_lock( Core::StateEngine::GetMutex() );	
	if(	PreferencesManager::Instance()->auto_save_state_->get() == false )
	{
		return false;
	}

	double timeout = PreferencesManager::Instance()->auto_save_time_state_->get() * 60;
	double time_duration = ProjectManager::Instance()->get_time_since_last_saved_session();

	if( time_duration < timeout )
	{
		return false;
	}
	
// 	if( PreferencesManager::Instance()->smart_save_state_->get() == true )
// 	{
// 		lock_type lock( this->get_mutex() );
// 		while( true )
// 		{
// 			if( !Core::ActionDispatcher::Instance()->is_busy() )
// 			{
// 				boost::posix_time::ptime last_action_completed = 
// 					Core::ActionDispatcher::Instance()->last_action_completed();
// 				boost::posix_time::ptime current_time = 
// 					boost::posix_time::second_clock::local_time();
// 				
// 				boost::posix_time::time_duration duration = current_time - last_action_completed;
// 				double time_since_last_action = static_cast< double >( duration.total_milliseconds() ) * 0.001;
// 
// 				if( time_since_last_action > 5  )
// 				{
// 					return true;
// 				}
// 			}
// 
// 			boost::posix_time::milliseconds wait_time( static_cast< int >( 5000.0 ) );
// 			recompute_auto_save_.timed_wait( lock, wait_time );
// 		}




/*	}*/

		//if( Core::ActionDispatcher::Instance()->get_time_since_last_action() < 5.0 )
		//{
		//	return false;
		//}
	

	return true;
}

void AutoSave::do_auto_save()
{
	ActionSaveSession::Dispatch( Core::Interface::GetWidgetActionContext(), true );
}


} // end namespace Seg3D