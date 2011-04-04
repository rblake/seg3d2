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


// Core Includes 
#include <Core/Utils/Log.h>
#include <Core/Math/MathFunctions.h>

// UI includes
#include "ui_QtSliderDoubleCombo.h"

// QtUtils includes
#include <QtUtils/Widgets/QtSliderDoubleCombo.h>

namespace QtUtils
{

class QtSliderDoubleComboPrivate
{
public:
    Ui::SliderDoubleCombo ui_;
	double min_;
	double max_;
};

QtSliderDoubleCombo::QtSliderDoubleCombo( QWidget* parent, bool edit_range ) :
	QWidget( parent ),
	value_( 0 ),
    private_( new QtSliderDoubleComboPrivate )
{
    this->private_->ui_.setupUi( this );

    this->connect( this->private_->ui_.horizontalSlider, SIGNAL( valueChanged( int )), 
		this, SLOT( slider_signal( int )) );
    this->connect( this->private_->ui_.spinBox, SIGNAL( valueChanged( double )), 
		this, SLOT( spinner_signal( double )) );

	this->private_->ui_.horizontalSlider->setPageStep( 10000 );
	
}

QtSliderDoubleCombo::~QtSliderDoubleCombo()
{
}
	
void QtSliderDoubleCombo::set_description( std::string description )
{
	this->private_->ui_.description_->setText( QString::fromStdString( description ) );
}

// signal from the spinner
void QtSliderDoubleCombo::spinner_signal( double value )
{   
	this->value_ = Core::Clamp( value, this->private_->min_, this->private_->max_ );

    this->private_->ui_.horizontalSlider->blockSignals( true );
    this->private_->ui_.horizontalSlider->setValue( Core::Round( this->value_ * 10000 ) );
    Q_EMIT valueAdjusted( this->value_ );
	this->private_->ui_.horizontalSlider->blockSignals( false );
}

// signal from the slider
void QtSliderDoubleCombo::slider_signal( int value )
{
    this->private_->ui_.spinBox->blockSignals( true );
	this->value_ = Core::Clamp( value / 10000.0, 
		this->private_->min_, this->private_->max_ );
    this->private_->ui_.spinBox->setValue( this->value_ );
    Q_EMIT valueAdjusted( this->value_ );
	this->private_->ui_.spinBox->blockSignals( false );
}

void QtSliderDoubleCombo::setStep( double step )
{
    this->block_signals( true );
    int int_step = static_cast< int >( step * 10000 );
    this->private_->ui_.horizontalSlider->setSingleStep( int_step );
    this->private_->ui_.spinBox->setSingleStep( step );
    this->block_signals( false );
}

void QtSliderDoubleCombo::setRange( double min, double max )
{
    this->block_signals( true );
	this->private_->min_ = min;
	this->private_->max_ = max;
    this->private_->ui_.horizontalSlider->setRange( static_cast<int>( min * 10000 ), 
		static_cast<int>( max * 10000 ) );
    this->private_->ui_.spinBox->setRange( min, max );
    this->private_->ui_.min_->setNum( min );
    this->private_->ui_.max_->setNum( max );
    
    double tick = ( ( max * 10000.0 ) - ( min * 10000.0 ) ) / 10.0;
    this->private_->ui_.horizontalSlider->setTickInterval( tick );
    this->block_signals( false );
}
void QtSliderDoubleCombo::setCurrentValue( double value )
{
    this->block_signals( true );
	this->value_ = Core::Clamp( value, this->private_->min_, this->private_->max_ );
    this->private_->ui_.horizontalSlider->setValue( static_cast<int>( this->value_ * 10000.0 ) );
    this->private_->ui_.spinBox->setValue( this->value_ );
    this->block_signals( false );
    Q_EMIT valueAdjusted( this->value_ );
}

void QtSliderDoubleCombo::change_min( double new_min )
{
	this->private_->min_ = new_min;
    this->block_signals( true );
    this->private_->ui_.horizontalSlider->setMinimum( static_cast<int>( new_min * 10000.0 ) );
    this->private_->ui_.spinBox->setMinimum( new_min );
    this->private_->ui_.min_->setNum( new_min );
    int tick = ( this->private_->ui_.max_->text().toInt() - 
		this->private_->ui_.min_->text().toInt()) / 10;
    this->private_->ui_.horizontalSlider->setTickInterval( tick * 10000 );
    this->block_signals( false );
}

void QtSliderDoubleCombo::change_max( double new_max )
{
	this->private_->max_ = new_max;
    this->block_signals( true );
    this->private_->ui_.horizontalSlider->setMaximum( static_cast<int>( new_max * 10000.0 ) );
    this->private_->ui_.spinBox->setMaximum( new_max );
    this->private_->ui_.max_->setNum( new_max );
    int tick = (this->private_->ui_.max_->text().toInt() - 
		this->private_->ui_.min_->text().toInt()) / 10;
    this->private_->ui_.horizontalSlider->setTickInterval( tick * 10000 );
    this->block_signals( false );
}

void QtSliderDoubleCombo::block_signals( bool block )
{
    this->private_->ui_.horizontalSlider->blockSignals( block );
    this->private_->ui_.spinBox->blockSignals( block ); 
}

}  // end namespace QtUtils
