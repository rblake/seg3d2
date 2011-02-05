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

#include <list>
#include <vector>

#include <boost/foreach.hpp>

#include <QCoreApplication>

#include <QtUtils/Widgets/QtTransferFunctionCurve.h>
#include <QtUtils/Widgets/QtTransferFunctionControlPoint.h>

namespace QtUtils
{

typedef std::list< QtTransferFunctionControlPoint* > control_point_list_type;
typedef std::vector< QGraphicsLineItem* > edge_list_type;

//////////////////////////////////////////////////////////////////////////
// Class QtTransferFunctionCurvePrivate
//////////////////////////////////////////////////////////////////////////

class QtTransferFunctionCurvePrivate
{
public:
	void update_edges();
	void update_curve_appearance();
	qreal get_edge_width();
	qreal get_point_radius();
	void get_control_points( Core::TransferFunctionControlPointVector& points );

	control_point_list_type control_points_;
	edge_list_type edges_;
	QGraphicsScene* scene_;
	bool active_;
	QColor color_;
	std::string feature_id_;
};

bool CompareControlPoints( const QtTransferFunctionControlPoint* const lhs,
	const QtTransferFunctionControlPoint* const rhs )
{
	return lhs->scenePos().x() < rhs->scenePos().x();
}

void QtTransferFunctionCurvePrivate::update_edges()
{
	if ( this->control_points_.size() > 1 )
	{
		this->control_points_.sort( CompareControlPoints );
		size_t edge_index = 0;
		control_point_list_type::const_iterator second_it = this->control_points_.begin();
		control_point_list_type::const_iterator first_it = second_it++;
		while ( second_it != this->control_points_.end() )
		{
			QGraphicsLineItem* edge = 0;
			if ( edge_index >= this->edges_.size() )
			{
				edge = new QGraphicsLineItem( 0, this->scene_ );
				this->edges_.push_back( edge );
				edge->setZValue( -1 );
			}
			else
			{
				edge = this->edges_[ edge_index ];
			}
			edge->setLine( QLineF( ( *first_it )->scenePos(), ( *second_it )->scenePos() ) );
			first_it = second_it++;
			++edge_index;
		}
	}

	while ( this->edges_.size() > 0 &&
		this->edges_.size() + 1 > this->control_points_.size() )
	{
		QGraphicsLineItem* extra_edge = this->edges_.back();
		this->scene_->removeItem( extra_edge );
		delete extra_edge;
		this->edges_.pop_back();
	}
}

qreal QtTransferFunctionCurvePrivate::get_edge_width()
{
	return this->active_ ? 3.0 : 2.0;
}

qreal QtTransferFunctionCurvePrivate::get_point_radius()
{
	return this->active_ ? 7.0 : 5.0;
}

void QtTransferFunctionCurvePrivate::update_curve_appearance()
{
	qreal edge_width = this->get_edge_width();
	qreal point_radius = this->get_point_radius();
	QPen edge_pen( this->color_, edge_width );
	edge_pen.setCosmetic( true );
	QBrush point_brush( this->color_ );

	control_point_list_type::const_iterator point_it = this->control_points_.begin();
	while ( point_it != this->control_points_.end() )
	{
		QtTransferFunctionControlPoint* point = *point_it;
		point->setRect( -point_radius, -point_radius, point_radius * 2 + 1, point_radius * 2 + 1 );
		point->setBrush( point_brush );
		++point_it;
	}

	edge_list_type::const_iterator edge_it = this->edges_.begin();
	while ( edge_it != this->edges_.end() )
	{
		( *edge_it )->setPen( edge_pen );
		++edge_it;
	}
}

void QtTransferFunctionCurvePrivate::get_control_points( 
	Core::TransferFunctionControlPointVector& points )
{
	QRectF scene_rect =  this->scene_->sceneRect();
	points.clear();
	points.reserve( this->control_points_.size() );
	BOOST_FOREACH( QtTransferFunctionControlPoint* control_point, this->control_points_ )
	{
		points.push_back( Core::TransferFunctionControlPoint( 
			( control_point->scenePos().x() - scene_rect.left() ) / scene_rect.width(),
			( scene_rect.bottom() - control_point->scenePos().y() ) / scene_rect.height() ) );
	}
}

//////////////////////////////////////////////////////////////////////////
// Class QtTransferFunctionCurve
//////////////////////////////////////////////////////////////////////////

QtTransferFunctionCurve::QtTransferFunctionCurve( const std::string& feature_id, 
												 QGraphicsScene* parent ) :
	QObject( parent ),
	private_( new QtTransferFunctionCurvePrivate )
{
	this->private_->scene_ = parent;
	this->private_->active_ = false;
	this->private_->color_ = QColor( 255, 255, 255 );
	this->private_->feature_id_ = feature_id;
}

QtTransferFunctionCurve::~QtTransferFunctionCurve()
{
	if ( !QCoreApplication::closingDown() )
	{
		BOOST_FOREACH( QtTransferFunctionControlPoint* control_point, this->private_->control_points_ )
		{
			this->private_->scene_->removeItem( control_point );
			delete control_point;
		}

		BOOST_FOREACH( QGraphicsLineItem* line_item, this->private_->edges_ )
		{
			this->private_->scene_->removeItem( line_item );
			delete line_item;
		}
	}

	delete this->private_;
}

void QtTransferFunctionCurve::add_control_point( const QPointF& pos )
{
	QtTransferFunctionControlPoint* point = new QtTransferFunctionControlPoint( this );
	point->setPos( pos );
	this->private_->scene_->addItem( point );

	this->connect( point, SIGNAL( position_changed( QPointF ) ),
		SLOT( handle_control_point_moved() ) );
	this->private_->control_points_.push_back( point );
	this->private_->update_edges();
	this->private_->update_curve_appearance();

	Core::TransferFunctionControlPointVector points;
	this->private_->get_control_points( points );
	Q_EMIT this->control_points_changed( points );
}

void QtTransferFunctionCurve::remove_control_point( QtTransferFunctionControlPoint* control_point )
{
	control_point->disconnect( this );
	this->private_->scene_->removeItem( control_point );
	this->private_->control_points_.remove( control_point );
	control_point->deleteLater();
	this->private_->update_edges();

	Core::TransferFunctionControlPointVector points;
	this->private_->get_control_points( points );
	Q_EMIT this->control_points_changed( points );
}

void QtTransferFunctionCurve::set_control_points( 
	const Core::TransferFunctionControlPointVector& points )
{
	size_t num_pts = points.size();

	while ( this->private_->control_points_.size() < num_pts )
	{
		QtTransferFunctionControlPoint* point_item = new QtTransferFunctionControlPoint( this );
		this->private_->scene_->addItem( point_item );
		this->connect( point_item, SIGNAL( position_changed( QPointF ) ),
			SLOT( handle_control_point_moved() ) );
		this->private_->control_points_.push_back( point_item );
	}

	while ( this->private_->control_points_.size() > num_pts )
	{
		QtTransferFunctionControlPoint* point_item = this->private_->control_points_.back();
		this->private_->control_points_.pop_back();
		this->private_->scene_->removeItem( point_item );
		delete point_item;
	}

	QRectF scene_rect = this->private_->scene_->sceneRect();
	control_point_list_type::iterator item_it = this->private_->control_points_.begin();
	for ( size_t i = 0; i < num_pts; ++i, ++item_it )
	{
		QPointF pt_pos( scene_rect.left() + scene_rect.width() * points[ i ].get_value(),
			scene_rect.bottom() - scene_rect.height() * points[ i ].get_opacity() );
		pt_pos.setX( qMin( scene_rect.right(), qMax( pt_pos.x(), scene_rect.left() ) ) );
		pt_pos.setY( qMin( scene_rect.bottom(), qMax( pt_pos.y(), scene_rect.top() ) ) );
		( *item_it )->setPos( pt_pos );
	}

	this->private_->update_edges();
	this->private_->update_curve_appearance();
}

void QtTransferFunctionCurve::handle_control_point_moved()
{
	this->private_->update_edges();

	Core::TransferFunctionControlPointVector points;
	this->private_->get_control_points( points );
	Q_EMIT this->control_points_changed( points );
}

void QtTransferFunctionCurve::set_active( bool active )
{
	if ( this->private_->active_ != active )
	{
		this->private_->active_ = active;
		this->private_->update_curve_appearance();
		Q_EMIT this->activated( active );
	}
}

void QtTransferFunctionCurve::set_color( const QColor& color )
{
	if ( this->private_->color_ != color )
	{
		this->private_->color_ = color;
		this->private_->update_curve_appearance();
	}
}

const std::string& QtTransferFunctionCurve::get_feature_id() const
{
	return this->private_->feature_id_;
}

bool QtTransferFunctionCurve::is_active() const
{
	return this->private_->active_;
}

} // end namespace QtUtils
