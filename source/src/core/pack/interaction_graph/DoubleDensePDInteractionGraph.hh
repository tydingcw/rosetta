// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/DoubleDensePDInteractionGraph.hh
/// @brief  Double-The-Memory interaction graph which stores interaction energies
/// for a rotamer and all of its neighbors in a single row in memory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pack_interaction_graph_DoubleDensePDInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_DoubleDensePDInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/DoubleDensePDInteractionGraph.fwd.hh>

// Package Headers

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>

//STL Headers

//ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

#include <ObjexxFCL/FArray2.hh>


namespace core {
namespace pack {
namespace interaction_graph {

class DoubleDensePDNode;
class DoubleDensePDEdge;
class DoubleDensePDInteractionGraph;

class DoubleDensePDNode : public PrecomputedPairEnergiesNode
{
public:
	DoubleDensePDNode(InteractionGraphBase * owner, int node_id, int num_states);
	~DoubleDensePDNode() override;
	void print() const override;

	void update_one_body_energy( int state, core::PackerEnergy energy ) override;
	virtual void update_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies );
	void add_to_one_body_energy( int state, core::PackerEnergy energy ) override;
	void add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energies ) override;
	void zero_one_body_energies() override;
	core::PackerEnergy get_one_body_energy( int state ) const override;

	void prepare_for_simulated_annealing() override;
	//virtual unsigned int getMemoryUsageInBytes() const;

	void assign_zero_state() override;
	bool state_unassigned() const override { return current_state_ == 0;}
	void assign_state(int new_state);
	int get_current_state() const;
	core::PackerEnergy get_one_body_energy_current_state();
	inline core::PackerEnergy project_deltaE_for_substitution
	(
		int alternate_state,
		core::PackerEnergy & prev_node_energy
	);
	void commit_considered_substitution();

	inline
	void acknowledge_neighbors_state_substitution(
		int edge_to_altered_neighbor,
		core::PackerEnergy new_edge_energy,
		int other_node_new_state);

	void print_internal_energies() const;

	void update_internal_energy_sums();

	unsigned int count_static_memory() const override;
	unsigned int count_dynamic_memory() const override;

protected:
	void update_internal_vectors();

private:

	inline
	DoubleDensePDEdge* get_incident_dpd_edge( int index ) const
	{
		return ( DoubleDensePDEdge* ) get_incident_edge( index ); // c-style cast since static_cast won't compile
	}

	inline
	DoubleDensePDNode* get_adjacent_dpd_node( int index ) const
	{
		return ( DoubleDensePDNode* ) get_adjacent_node( index );
	}

	inline
	DoubleDensePDInteractionGraph* get_dpdig_owner() const
	{
		return ( DoubleDensePDInteractionGraph* ) get_owner();  // c-style cast since static_cast won't compile
	}

	std::vector< core::PackerEnergy > one_body_energies_;

	std::vector< int > neighbors_curr_state_;
	std::vector< int > neighbors_curr_state_plus_offset_;
	std::vector< int > neighbors_num_states_;
	std::vector< int > neighbors_rotindex_offset_;
	ObjexxFCL::FArray2D< core::PackerEnergy > rotamer_energies_; // dim1: index of neighbor state, dim2 index of my state

	int current_state_;
	core::PackerEnergy curr_state_one_body_energy_;
	core::PackerEnergy curr_state_total_energy_;
	std::vector< core::PackerEnergy > curr_state_two_body_energies_;

	int alternate_state_;
	core::PackerEnergy alternate_state_one_body_energy_;
	core::PackerEnergy alternate_state_total_energy_;
	std::vector< core::PackerEnergy > alternate_state_two_body_energies_;

	bool alternate_state_is_being_considered_;

	//no default constructor, uncopyable
	DoubleDensePDNode();
	DoubleDensePDNode( DoubleDensePDNode const & );
	DoubleDensePDNode & operator = ( DoubleDensePDNode const & );
};

class DoubleDensePDEdge : public PrecomputedPairEnergiesEdge
{
public:
	DoubleDensePDEdge(InteractionGraphBase* owner, int first_node_ind, int second_node_ind);
	~DoubleDensePDEdge() override;
	void set_sparse_aa_info(ObjexxFCL::FArray2_bool const & ) override {}
	bool get_sparse_aa_info( int, int) const override {return true;} //"all amino acids are neighbors"
	void add_to_two_body_energy(int const, int const, core::PackerEnergy const) override;
	void add_to_two_body_energies( ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array ) override;
	void set_two_body_energy(int const, int const, core::PackerEnergy const) override;
	void clear_two_body_energy(int const, int const) override;
	core::PackerEnergy get_two_body_energy( int const, int const ) const override;

	void force_aa_neighbors(int, int) override {} //all aa's are already neighbors -- dense representation
	void force_all_aa_neighbors() override {} //same thing


	void declare_energies_final() override;
	void prepare_for_simulated_annealing() override;
	//virtual unsigned int getMemoryUsageInBytes() const;

	core::PackerEnergy get_current_two_body_energy();

	void acknowledge_state_change(
		int node_ind,
		int new_state,
		core::PackerEnergy & new_energy
	);
	void acknowledge_state_zeroed( int node_ind );

	static
	inline
	core::PackerEnergy get_alternate_state_energy(
		int first_node_state,
		int second_node_state,
		ObjexxFCL::FArray2< core::PackerEnergy > & edge_energy_table
	);

	inline void acknowledge_substitution(
		int substituted_node_index,
		core::PackerEnergy const curr_state_energy,
		int nodes_new_state
	);

	int get_two_body_table_size() const;
	//ObjexxFCL::FArray2Da< core::PackerEnergy > get_edge_table_ptr();

	unsigned int count_static_memory() const override;
	unsigned int count_dynamic_memory() const override;

	void set_edge_weight( Real weight ) override;

private:
	inline
	DoubleDensePDNode* get_dpd_node( int index ) const
	{
		return (DoubleDensePDNode*) get_node( index ); //c-style cast since static_cast won't compile
	}

	inline
	DoubleDensePDInteractionGraph* get_dpdig_owner() const
	{
		return (DoubleDensePDInteractionGraph*) get_owner(); //c-style cast since static_cast won't compile
	}

	ObjexxFCL::FArray2D< core::PackerEnergy > two_body_energies_; //Dense matrix
	core::PackerEnergy curr_state_energy_;
	bool energies_updated_since_last_prep_for_simA_;

	//no default constructor, uncopyable
	DoubleDensePDEdge();
	DoubleDensePDEdge( DoubleDensePDEdge const & );
	DoubleDensePDEdge & operator = ( DoubleDensePDEdge const & );
};

class DoubleDensePDInteractionGraph : public PrecomputedPairEnergiesInteractionGraph
{
public:
	DoubleDensePDInteractionGraph(int num_nodes);
	void initialize( pack_basic::RotamerSetsBase const & rot_sets ) override;

	core::PackerEnergy get_one_body_energy_for_node_state( int node, int state) override;

	//virtual void set_num_aatypes(int) {}
	int  get_num_aatypes() const override {return 1;}

	void blanket_assign_state_0() override;
	core::PackerEnergy set_state_for_node(int node_ind, int new_state) override;
	core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states) override;
	void consider_substitution
	(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	) override;

	/// @brief Accepts (commits) the state change previously considered in a call to
	/// consider_substitution and returns the energy of the entire graph
	core::PackerEnergy commit_considered_substitution() override;

	/// @brief removes all accumulated numerical drift and returns the
	/// energy for the current state assignment.
	core::PackerEnergy get_energy_current_state_assignment() override;
	/// @brief returns the number of floats used in all edge two-body energy tables
	int get_edge_memory_usage() const override;

	/// @brief outputs the current state for each node, useful for debugging
	void print_current_state_assignment() const override;
	void set_errorfull_deltaE_threshold( core::PackerEnergy ) override {};

	/// @brief a user may define subsets of the vertex set for which they would like to
	/// know the internal energy sum.
	core::PackerEnergy get_energy_sum_for_vertex_group( int group_id ) override;

	unsigned int count_static_memory() const override;
	unsigned int count_dynamic_memory() const override;

protected:
	//virtual unsigned int getMemoryUsageInBytes() const;

	NodeBase* create_new_node( int node_index, int num_states) override;
	EdgeBase* create_new_edge( int index1, int index2) override;

	/// @brief removes numerical drift that can accumulate over the course of
	/// many state assignment changes within simulated annealing
	void update_internal_energy_totals();

	inline
	DoubleDensePDNode* get_dpd_node(int index) const
	{
		return (DoubleDensePDNode*) get_node( index );
	}

private:
	int num_commits_since_last_update_;
	core::PackerEnergy total_energy_current_state_assignment_;
	core::PackerEnergy total_energy_alternate_state_assignment_;
	int node_considering_alt_state_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10

	//no default constructor, uncopyable
	DoubleDensePDInteractionGraph();
	DoubleDensePDInteractionGraph( DoubleDensePDInteractionGraph const & );
	DoubleDensePDInteractionGraph & operator = ( DoubleDensePDInteractionGraph const & );
};

/// @brief returns the change in energy that would be induced by switching this node
/// from its current state into another state
///
/// iterates across the incident edges for a node in two phases:
/// in the first phase, it examines edges leading to higher-indexed nodes
/// in the second phase, it examines edges leading to smaller-indexed nodes.
/// for cache efficiency, all of the amino-acid-neighbor-offset information
/// that each edge calculates is stored on the nodes themselves.  The edges
/// are never touched; rather, their private information is stored on the nodes
/// and handed to static member functions of the DoubleDensePDEdge class.  This "store
/// edge information on the nodes" strategy gives me performance equivalent
/// to the previous energy2b lookup tables.
///
/// @param alternate_state - [in] - the alternate state to consider
/// @param previous_energy_for_node - [out] - the old energy1b/energy2b sum for this
/// node; used by simulate annealing.
///
inline
core::PackerEnergy
DoubleDensePDNode::project_deltaE_for_substitution(
	int alternate_state,
	core::PackerEnergy & prev_energy_for_node
)
{

	alternate_state_is_being_considered_ = true;
	//std::cerr << "proj_deltaE: node -  " << get_node_index()
	//<< " alt state " << alternate_state << "...";

	alternate_state_ = alternate_state;
	alternate_state_one_body_energy_ = one_body_energies_[ alternate_state ];
	alternate_state_total_energy_ = alternate_state_one_body_energy_;
	prev_energy_for_node = curr_state_total_energy_;

	if ( get_num_incident_edges() == 0 ) {
		return alternate_state_total_energy_ - curr_state_total_energy_;
	}
	int const altstate_offset = rotamer_energies_.index( 1, alternate_state_ ) - 1;

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		alternate_state_two_body_energies_[ ii ] = rotamer_energies_[
			altstate_offset + neighbors_curr_state_plus_offset_[ ii ] ];
		//alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
		//std::cerr << " edge " << ii << " E= " << alternate_state_two_body_energies_[ ii ];
	}

	for ( int ii = 1; ii <= get_num_incident_edges(); ++ii ) {
		//alternate_state_two_body_energies_[ ii ] = rotamer_energies_[
		// altstate_offset + neighbors_curr_state_plus_offset_[ ii ] ];
		alternate_state_total_energy_ += alternate_state_two_body_energies_[ ii ];
		//std::cerr << " edge " << ii << " E= " << alternate_state_two_body_energies_[ ii ];
	}

	//std::cerr<< "..done" << std::endl;

	return alternate_state_total_energy_ - curr_state_total_energy_;

}

/// @brief updates bookkeeping arrays for when a neighbor has changed its state
///
/// @param edge_to_altered_neighbor - [in] - the index for the edge that connects
///  this node to the node that just changed its state
/// @param new_edge_energ - [in] - the pair energy between this node in its current
/// state and the new state of the node that just changed its state
/// @param other_node_new_state - [in] - the state the neighbor just adopted
inline
void DoubleDensePDNode::acknowledge_neighbors_state_substitution(
	int edge_to_altered_neighbor,
	core::PackerEnergy new_edge_energy,
	int other_node_new_state
)
{

	curr_state_total_energy_ +=
		new_edge_energy - curr_state_two_body_energies_[edge_to_altered_neighbor];
	curr_state_two_body_energies_[edge_to_altered_neighbor] = new_edge_energy;
	neighbors_curr_state_[ edge_to_altered_neighbor ] = other_node_new_state;
	neighbors_curr_state_plus_offset_[ edge_to_altered_neighbor ] = other_node_new_state + neighbors_rotindex_offset_[ edge_to_altered_neighbor ];
	return;
}

/// @brief static method that looks up the two body energy when the
/// node with the smaller index on an edge is considering an alternate state
///
/// @param first_node_alt_state - [in] - the alternate state for the lower-indexed node
/// @param second_node_orig_state - [in] - the current state for the higher-indexed node
/// @param edge_energy_table - [in] - the proxy FArray pointing at the edge table
///  connecting the two nodes.
inline
core::PackerEnergy
DoubleDensePDEdge::get_alternate_state_energy(
	int first_node_state,
	int second_node_state,
	ObjexxFCL::FArray2< core::PackerEnergy > & edge_energy_table
)
{
	if ( first_node_state == 0 || second_node_state == 0 ) {
		return 0.0f;
	} else {
		return edge_energy_table( second_node_state, first_node_state );
	}
}


/// @brief update bookkeeping information when one of the nodes an edge is incident
/// upon changes state
///
/// @param substituted_node_index - [in] - index of the node that chagned its state
/// @param curr_state_energy - [in] - the two body energy given the new state
/// @param nodes_new_state - [in] - the state the node just transitioned into
/// @param nodes_new_state_sparse_info - [in] - sparse matrix info for the new state
inline
void
DoubleDensePDEdge::acknowledge_substitution(
	int substituted_node_index,
	core::PackerEnergy const curr_state_energy,
	int nodes_new_state
)
{
	int node_substituted = substituted_node_index == get_node_index(0) ? 0 : 1;
	int node_not_substituted = ! node_substituted;

	curr_state_energy_ = curr_state_energy;

	get_dpd_node( node_not_substituted )->
		acknowledge_neighbors_state_substitution
		(
		get_edges_position_in_nodes_edge_vector( node_not_substituted ),
		curr_state_energy_,
		nodes_new_state
	);
}


} // namespace interaction_graph
} // namespace pack
} // namespace core

#endif
