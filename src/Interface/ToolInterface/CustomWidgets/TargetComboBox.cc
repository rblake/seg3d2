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

// Interface includes
#include <Interface/ToolInterface/CustomWidgets/TargetComboBox.h>

// Application Includes
#include <Application/Layer/Layer.h>
#include <Application/LayerManager/LayerManager.h>

namespace Seg3D
{
	
	
TargetComboBox::TargetComboBox( QWidget *parent )
{
	this->setParent( parent );
	this->setMinimumHeight( 26 );
	
	QSizePolicy size_policy(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
	this->setSizePolicy( size_policy );
	
	this->add_connection( LayerManager::Instance()->layers_changed_signal_.connect(
		boost::bind( &TargetComboBox::sync_layers, this ) ) );
	
	this->sync_layers();
}

TargetComboBox::~TargetComboBox()
{
	this->disconnect_all();
}
	
void TargetComboBox::sync_layers()
{
	this->value_ = this->currentText().toStdString();
	std::vector< LayerHandle > target_layers;
	LayerManager::Instance()->get_layers( target_layers );
	bool has_layers = false;
	
	//this->value_ = this->currentText().toStdString();
		
	this->clear();
	for( int i = ( static_cast< int >( target_layers.size() ) - 1 ); i > -1; i-- )
	{   
		has_layers = true;
		if ( target_layers[i]->type() == Core::VolumeType::DATA_E ) 
		{
			this->addItem( QString::fromStdString( target_layers[i]->get_layer_name() ) );
		}	
	} 
	
	if( this->value_ != "" ) 
	{
		int index = this->findText( QString::fromStdString( this->value_ ), 
			Qt::MatchFlags( Qt::CaseInsensitive ) );
		this->setCurrentIndex( index );
	}
	
	Q_EMIT valid( has_layers );
	
	this->update();
}	

} // end namespace Seg3D