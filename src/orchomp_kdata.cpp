/** \file orchomp_kdata.cpp
 * \brief Implementation of orchomp_rdata, a parser for sphere data
 *        provided with an OpenRAVE kinbody XML file.
 * \author Christopher Dellin
 * \date 2012
 */

/* (C) Copyright 2012-2013 Carnegie Mellon University */

/* This module (orchomp) is part of libcd.
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

#include <openrave/openrave.h>

#include "orchomp_kdata.h"

namespace orchomp
{

kdata::kdata() : OpenRAVE::XMLReadable("orchomp")
{
   this->sphereelems = 0;
}

kdata::~kdata()
{
   struct sphereelem * e;
   while (this->sphereelems)
   {
      e=this->sphereelems->next;
      free(this->sphereelems->s);
      free(this->sphereelems);
      this->sphereelems = e;
   }
}


kdata_parser::kdata_parser(boost::shared_ptr<kdata> passed_d, const OpenRAVE::AttributesList& atts)
{
   /* save or construct the kdata object */
   this->d = passed_d;
   if(!this->d) this->d.reset(new kdata());
   /* get ready */
   this->inside_spheres = false;
}

OpenRAVE::XMLReadablePtr kdata_parser::GetReadable()
{
   return this->d;
}

OpenRAVE::BaseXMLReader::ProcessElement kdata_parser::startElement(const std::string& name, const OpenRAVE::AttributesList& atts)
{
   if (name == "spheres")
   {
      if (this->inside_spheres) RAVELOG_ERROR("you can't have <spheres> inside <spheres>!\n");
      this->inside_spheres = true;
      return PE_Support;
   }
   if (name == "sphere")
   {
      struct sphere * s;
      struct sphereelem * e;
      if (!this->inside_spheres) { RAVELOG_ERROR("you can't have <sphere> not inside <spheres>!\n"); return PE_Pass; }
      s = (struct sphere *) malloc(sizeof(struct sphere));
      for(OpenRAVE::AttributesList::const_iterator itatt = atts.begin(); itatt != atts.end(); ++itatt)
      {
         if (itatt->first=="link")
            strcpy(s->linkname, itatt->second.c_str());
         else if (itatt->first=="radius")
            s->radius = strtod(itatt->second.c_str(), 0);
         else if (itatt->first=="pos")
            sscanf(itatt->second.c_str(), "%lf %lf %lf", &s->pos[0], &s->pos[1], &s->pos[2]);
         else
            RAVELOG_ERROR("unknown attribute %s=%s!\n",itatt->first.c_str(),itatt->second.c_str());
      }
      /* insert at head of kdata list */
      e = (struct sphereelem *) malloc(sizeof(struct sphereelem));
      e->s = s;
      e->next = this->d->sphereelems;
      this->d->sphereelems = e;
      return PE_Support;
   }
   return PE_Pass;
}

void kdata_parser::characters(const std::string& ch)
{
   return;
}

bool kdata_parser::endElement(const std::string& name)
{
   if (name == "orchomp") return true;
   if (name == "spheres")
   {
      if (!this->inside_spheres) RAVELOG_ERROR("you can't have </spheres> without matching <spheres>!\n");
      this->inside_spheres = false;
   }
   else if (name == "sphere")
   {
      if (!this->inside_spheres) RAVELOG_ERROR("you can't have </sphere> not inside <spheres>!\n");
   }
   else
      RAVELOG_ERROR("unknown field %s\n", name.c_str());
   return false;
}

} /* namespace orchomp */
