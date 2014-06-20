/** \file orcdchomp_kdata.h
 * \brief Interface to orcdchomp_kdata, a parser for sphere data provided
 *        with an OpenRAVE kinbody XML file.
 * \author Christopher Dellin
 * \date 2012
 */

/* (C) Copyright 2012-2013 Carnegie Mellon University */

/* This module (orcdchomp) is part of libcd.
 *
 * This module of libcd is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This module of libcd is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public License is provided with libcd
 * (license-gpl.txt) and is also available at <http://www.gnu.org/licenses/>.
 */

/* requires:
 *  - openrave/openrave.h
 * */

#ifndef _ORCHOMP_KDATA_H_
#define _ORCHOMP_KDATA_H_


namespace orchomp
{


class Sphere
{
  public:
    // The radius of the sphere
    double radius;

    //a pointer to the kinBody that
    //the sphere comes off of.
    OpenRAVE::KinBody::Link * link;
    
    OpenRAVE::KinBody * body;

    //the name of the kinbody it is linked to;
    std::string linkname; 
    
    //the index of the robot link
    int linkindex;

    //The transform between the coordinates of the
    //  robot link and the sphere. Since it is a sphere,
    //  rotation is meaningless, so it should just be an
    //  xyz vector.
    double pose[3];

    Sphere() : link( NULL ), linkindex(-1){}
};


/* the kinbody-attached data class */
class kdata : public OpenRAVE::XMLReadable
{
public:
   std::vector<Sphere> spheres;
   kdata();
   ~kdata();
};


/* the kdata-parser */
class kdata_parser : public OpenRAVE::BaseXMLReader
{
public:
   boost::shared_ptr<kdata> d;
   bool inside_spheres;

   kdata_parser(boost::shared_ptr<kdata> passed_d, const OpenRAVE::AttributesList& atts);
   virtual OpenRAVE::XMLReadablePtr GetReadable();
   virtual ProcessElement startElement(const std::string& name, const OpenRAVE::AttributesList& atts);
   virtual void characters(const std::string& ch);
   virtual bool endElement(const std::string& name);
};

} /* namespace orcdchomp */



#endif
