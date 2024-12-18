// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/FilterFactory.cc
/// @brief
/// @author ashworth

#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/filter_schemas.hh>

// Package headers
#include <protocols/filters/BasicFilters.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Basic headers
#include <basic/citation_manager/CitationManager.hh>

// Boost headers

namespace protocols {
namespace filters {


static basic::Tracer TR( "protocols.filters.FilterFactory" );

FilterFactory::FilterFactory() = default;

FilterFactory::~FilterFactory()= default;

/// @brief add a Filter prototype, using its default type name as the map key
void
FilterFactory::factory_register( FilterCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const filter_type( creator->keyname() );
	if ( filter_type == "UNDEFINED NAME" ) {
		utility_exit_with_message("Can't map derived Filter with undefined type name.");
	}
	if ( filter_creator_map_.find( filter_type ) != filter_creator_map_.end() ) {
		utility_exit_with_message("FilterFactory::factory_register already has a filter creator with name \"" + filter_type + "\".  Conflicting Filter names" );
	}
	filter_creator_map_[ filter_type ] = creator;
}

/// @brief Is there a filter with the given name that's known to Rosetta?
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
FilterFactory::filter_exists(
	std::string const &filter_name
) const {
	return( filter_creator_map_.find( filter_name ) != filter_creator_map_.end() );
}

/// @brief Get the XML schema for a given filter.
/// @details Throws an error if the filter is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
FilterFactory::provide_xml_schema(
	std::string const &filter_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	auto iter( filter_creator_map_.find( filter_name ) );
	if ( iter != filter_creator_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: FilterCreatorOP prototype for " + filter_name + " is NULL!" );
		}
		iter->second->provide_xml_schema( xsd );
	} else {
		TR << "Available filters: ";
		for ( auto const & filt_it : filter_creator_map_ ) {
			TR << filt_it.first << ", ";
		}
		TR << std::endl;
		utility_exit_with_message( filter_name + " is not known to the FilterFactory. Was it registered in the appropriate initialization files (src/protocols/init/init.FilterCreators.ihh and src/protocols/init/init.FilterRegistrators.ihh)??" );
	}
}


/// @brief return new Filter by key lookup in filter_prototype_map_ (new Filter parses Tag if provided)
FilterOP
FilterFactory::newFilter( std::string const & filter_type )
const
{
	FilterMap::const_iterator iter( filter_creator_map_.find( filter_type ) );
	if ( iter != filter_creator_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: FilterCreatorOP prototype for " + filter_type + " is NULL!" );
		}
		// use of cloning method would be faithful to pre-initialized prototypes
		//return iter->second->clone();
		// fresh_instance prevents propagation of pre-initialized prototypes, which may be safer(?)
		return iter->second->create_filter();
	} else {
		TR<<"Available filters: ";
		for ( FilterMap::const_iterator filt_it = filter_creator_map_.begin(); filt_it != filter_creator_map_.end(); ++filt_it ) {
			TR<<filt_it->first<<", ";
		}
		TR<<std::endl;
		utility_exit_with_message( filter_type + " is not known to the FilterFactory. Was it registered in the appropriate initialization files (src/protocols/init/init.FilterCreators.ihh and src/protocols/init/init.FilterRegistrators.ihh)?" );
		return nullptr;
	}
}

/// @brief return new Filter by Tag parsing
/*FilterOP
FilterFactory::newFilter(
TagCOP const tag,
basic::datacache::DataMap & data
) {
FilterOP filter( newFilter( tag->getName() ) );
runtime_assert( filter );
filter->parse_my_tag( tag, data, filters, movers, pose );
return filter;
}*/

/// @brief return new Filter by Tag parsing
FilterOP
FilterFactory::newFilter(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	std::string user_defined_name
) const {
	FilterOP filter( newFilter( tag->getName() ) );
	runtime_assert( filter != nullptr );
	std::string filter_name;
	if ( tag->hasOption("name") ) {
		filter_name = tag->getOption<std::string>("name");
	} else if ( user_defined_name.size() != 0 ) {
		filter_name = user_defined_name;
	} else {
		utility_exit_with_message("Can't define unnamed Filter");
	}
	filter->set_user_defined_name( filter_name );
	filter->parse_my_tag( tag, data );

	//Register with the citation manager:
	basic::citation_manager::CitationCollectionList citations;
	filter->provide_citation_info( citations );
	basic::citation_manager::CitationManager::get_instance()->add_citations( citations );

	// if confidence specified, link to StochasticFilter and wrap inside CompoundFilter
	core::Real const confidence( tag->getOption< core::Real >( "confidence", 1.0 ) );
	if ( confidence <= 0.999 ) { // fuzzy logic
		// The split is `(1-confidence)` always true, `confidence` run subfilter.
		// Which means that StochasticFilter should be in a true state `(1-confidence)` of the time,
		// and run the subfilter when it's false.
		FilterOP stochastic_filter( new StochasticFilter( (1.0-confidence), filter->clone(), /*run_subfilter_on*/ false ) );
		stochastic_filter->set_user_defined_name( tag->getOption<std::string>("name") );

		basic::citation_manager::CitationCollectionList citations;
		stochastic_filter->provide_citation_info( citations );
		basic::citation_manager::CitationManager::get_instance()->add_citations( citations );

		return stochastic_filter;
	}
	return filter;
}


FilterFactory::FilterMap const & FilterFactory::filter_creator_map() const { return filter_creator_map_; }

void FilterFactory::define_filter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			filter_creator_map_,
			filter_xml_schema_group_name(),
			& complex_type_name_for_filter,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for Filters from FilterFactory; offending class"
			" must call protocols::filters::complex_type_name_for_filter when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::string FilterFactory::filter_xml_schema_group_name()
{
	return "filter";
}

/// @brief Get a human-readable listing of the citations for a given filter, by filter name.
/// @details Returns an empty string if there are no citations.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
std::string
FilterFactory::get_citation_humanreadable(
	std::string const & filter_name
) const {
	using namespace basic::citation_manager;
	CitationCollectionList citations;
	FilterOP filter( newFilter( filter_name ) );
	runtime_assert_string_msg( filter != nullptr, "Error in FilterFactory::get_citation_humanreadable(): Could not instantiate " + filter_name + "!" );
	filter->provide_citation_info(citations);
	if ( citations.empty() ) return "";
	std::ostringstream ss;
	ss << "References and author information for the " << filter_name << " filter:" << std::endl;
	ss << std::endl;
	basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info_from_list_to_stream( citations, ss );
	return ss.str();
}

} //namespace filters
} //namespace protocols
