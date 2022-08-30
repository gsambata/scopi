#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <vector>

#include <xtensor/xadapt.hpp>

#include "objects/types/base.hpp"
#include "types.hpp"
#include "property.hpp"
#include "crtp.hpp"

namespace scopi
{

    ////////////////////////////////
    // scopi_container definition //
    ////////////////////////////////
    /**
     * @brief Main data structure used in SCoPI.
     *
     * Array of particles.
     *
     * Inactive particles are placed at the begining of the container.
     *
     * In the following, "particle" means a base object (sphere, superellipsoid or plan) and an "object" can be a more complex object, such as a worm.
     * Particles are objects.
     *
     * @tparam dim Dimension (2 or 3).
     */
    template<std::size_t dim>
    class scopi_container
    {
    public:

        /**
         * @brief Alias for the type of the position.
         */
        using position_type = type::position_t<dim>;
        /**
         * @brief Alias for the type of the velocity.
         */
        using velocity_type = type::velocity_t<dim>;
        /**
         * @brief Alias for the type of the rotation vector (scalar in 2D, vector in 3D).
         */
        using rotation_type = type::rotation_t<dim>;
        /**
         * @brief Alias for the type of the force.
         */
        using force_type = type::force_t<dim>;
        /**
         * @brief Alias for the type of the mass.
         */
        using mass_type = double;
        /**
         * @brief Alias for the type of the momentum of inertia.
         */
        using moment_type = type::moment_t<dim>;
        /**
         * @brief Alias for the type of the quaternion.
         */
        using quaternion_type = type::quaternion_t;

        /**
         * @brief Constructor.
         */
        scopi_container();

        /**
         * @brief Reconstructs an object.
         *
         * An object can be a sphere, a superellipsoid, a plan, a worm,...
         *
         * @param i Index of the object.
         *
         * @return Object.
         */
        std::unique_ptr<object<dim, false>> operator[](std::size_t i);

        /**
         * @brief Appends the given element value to the end of the container. 
         *
         * @param s [in] Object to append. 
         * @param p [in] Properties of the object (see property.hpp).
         */
        void push_back(const object<dim>& s, const property<dim>& p = property<dim>());
        /**
         * @brief Copy a particle that was already in the container with a new position.
         *
         * The original particle remains unchanged.
         * All the other physical properties of the particle are copied, only the position is different.
         * Used for periodic boundary conditions.
         *
         * @param i [in] Index of the particle to copy.
         * @param pos [in] Position of the new particle.
         */
        void push_back(std::size_t i, const position_type& pos);

        /**
         * @brief Increase the capacity of the container.
         *
         * Increase the capacity of the container (the total number of particles that the container can hold without requiring reallocation) to a value that's greater or equal to \c size.
         * If \c size is greater than the current capacity, new storage is allocated, otherwise the function does nothing.
         *
         * \c reserve does not change the size of the container.
         *
         * @param size [in] New capacity of the container.
         */
        void reserve(std::size_t size);

        /**
         * @brief Array of particles' positions.
         */
        auto pos() const;
        /**
         * @brief Array of particles' positions.
         */
        auto pos();

        /**
         * @brief Array of particles' quaternions.
         */
        auto q() const;
        /**
         * @brief Array of particles' quaternions.
         */
        auto q();

        /**
         * @brief Array of particles' forces.
         */
        auto f() const;
        /**
         * @brief Array of particles' forces.
         */
        auto f();

        /**
         * @brief Array of particles' masses.
         */
        auto m() const;
        /**
         * @brief Array of particles' masses.
         */
        auto m();

        /**
         * @brief Array of particles' moments of inertia.
         */
        auto j() const;
        /**
         * @brief Array of particles' moments of inertia.
         */
        auto j();

        /**
         * @brief Array of particles' velocities.
         */
        auto v() const;
        /**
         * @brief Array of particles' velocities.
         */
        auto v();

        /**
         * @brief Array of particles' rotation.
         */
        auto omega() const;
        /**
         * @brief Array of particles' rotation.
         */
        auto omega();

        /**
         * @brief Array of particles' desired rotation.
         */
        auto desired_omega() const;
        /**
         * @brief Array of particles' desired rotations.
         */
        auto desired_omega();

        /**
         * @brief Array of particles' desired velocities.
         */
        auto vd() const;
        /**
         * @brief Array of particles' desired velocities.
         */
        auto vd();

        /**
         * @brief Number of objects in the container.
         */
        std::size_t size() const;
        /**
         * @brief Number of active particles in the container.
         */
        std::size_t nb_active() const;
        /**
         * @brief Number of inactive particles in the container.
         */
        std::size_t nb_inactive() const;

        /**
         * @brief Convert a particle index into an object index.
         *
         * The particle \c i is part of the object \c j, where \c j is the index of an object.
         * TODO a diagram would help a lot, but I didn't manage to do it.
         *
         * @param i [in] Index of a particle.
         *
         * @return Index of the object.
         */
        std::size_t object_index(std::size_t i) const;
        /**
         * @brief Convert an object index into a particle index.
         *
         * @param i [in] Index of an object.
         *
         * @return Index of the first particle in the object.
         */
        std::size_t offset(std::size_t i) const;

        /**
         * @brief Remove all fictive particles.
         */
        void reset_periodic();

    private:

        std::map<std::size_t, std::unique_ptr<base_constructor<dim>>> m_shape_map;
        std::vector<position_type> m_positions;  // pos()
        std::vector<quaternion_type> m_quaternions;  // q()
        std::vector<force_type> m_forces;  // f()
        std::vector<mass_type> m_masses;  // f()
        std::vector<moment_type> m_moments_inertia;  // f()
        std::vector<velocity_type> m_velocities;  // v()
        std::vector<velocity_type> m_desired_velocities;  // vd()
        std::vector<rotation_type> m_omega;  // omega()
        std::vector<rotation_type> m_desired_omega;  // desired_omega()
        std::vector<std::size_t> m_shapes_id;
        std::vector<std::size_t> m_offset;
        std::vector<std::size_t> m_periodic_indices;

        std::size_t m_periodic_ptr;

        std::size_t m_nb_inactive_core_objects;

        bool m_periodic_added;
    };

    template<std::size_t dim>
    scopi_container<dim>::scopi_container()
    : m_periodic_ptr(0)
    , m_nb_inactive_core_objects(0)
    , m_periodic_added(false)
    {}

    template<std::size_t dim>
    std::unique_ptr<object<dim, false>> scopi_container<dim>::operator[](std::size_t i)
    {
        return (*m_shape_map[m_shapes_id[i]])(&m_positions[m_offset[i]], &m_quaternions[m_offset[i]]);
    }

    template<std::size_t dim>
    void scopi_container<dim>::push_back(const object<dim>& s, const property<dim>& p)
    {
        assert(!m_periodic_added);

        if (m_offset.empty())
        {
            m_offset = {0, s.size()};
        }
        else
        {
            m_offset.push_back(m_offset.back() + s.size());
        }

        for(std::size_t i = 0; i< s.size(); ++i)
        {
            m_positions.push_back(s.pos(i));
            m_quaternions.push_back(s.q(i));
            m_velocities.push_back(p.velocity());
            m_omega.push_back(p.omega());
            m_desired_omega.push_back(p.desired_omega());
            m_desired_velocities.push_back(p.desired_velocity());
            m_forces.push_back(p.force());
            m_masses.push_back(p.mass());
            m_moments_inertia.push_back(p.moment_inertia());
        }

        if (!p.is_active())
        {
            m_nb_inactive_core_objects += s.size();
            if (m_nb_inactive_core_objects != m_positions.size())
            {
                throw std::runtime_error("All the obstacles must be pushed before the active particles.");
            }
        }

        auto it = m_shape_map.find(s.hash());
        if (it == m_shape_map.end())
        {
            m_shape_map.insert(std::make_pair(s.hash(), std::move(s.construct())));
        }

        m_shapes_id.push_back(s.hash());
        m_periodic_ptr++;
    }

    template<std::size_t dim>
    void scopi_container<dim>::push_back(std::size_t i,
                                         const position_type& pos)
    {
        assert(i >= 0 && i < m_positions.size());
        m_periodic_added = true;
        m_offset.push_back(m_offset.back() + m_offset[i+1] - m_offset[i]);

        m_positions.push_back(pos);
        m_quaternions.push_back(m_quaternions[i]);
        m_velocities.push_back(m_velocities[i]);
        m_omega.push_back(m_omega[i]);
        m_desired_omega.push_back(m_desired_omega[i]);
        m_desired_velocities.push_back(m_desired_velocities[i]);
        m_forces.push_back(m_forces[i]);
        m_masses.push_back(m_masses[i]);
        m_moments_inertia.push_back(m_moments_inertia[i]);

        m_shapes_id.push_back(m_shapes_id[i]);

        m_periodic_indices.push_back(i);
    }

    template<std::size_t dim>
    void scopi_container<dim>::reserve(std::size_t size)
    {
        m_positions.reserve(size);
        m_quaternions.reserve(size);
        m_velocities.reserve(size);
        m_desired_velocities.reserve(size);
        m_omega.reserve(size);
        m_desired_omega.reserve(size);
        m_forces.reserve(size);
        m_masses.reserve(size);
        m_moments_inertia.reserve(size);
        m_offset.reserve(size+1);
        m_shapes_id.reserve(size);
    }

    template<std::size_t dim>
    std::size_t scopi_container<dim>::size() const
    {
        return m_shapes_id.size();
    }

    template<std::size_t dim>
    std::size_t scopi_container<dim>::nb_active() const
    {
        return m_positions.size() - m_nb_inactive_core_objects;
    }

    template<std::size_t dim>
    std::size_t scopi_container<dim>::nb_inactive() const
    {
        return m_nb_inactive_core_objects;
    }

    // position

    template<std::size_t dim>
    auto scopi_container<dim>::pos() const
    {
        return xt::adapt(reinterpret_cast<const position_type*>(m_positions.data()), {m_positions.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::pos()
    {
        return xt::adapt(reinterpret_cast<position_type*>(m_positions.data()), {m_positions.size()});
    }

    // rotation

    template<std::size_t dim>
    auto scopi_container<dim>::q() const
    {
        return xt::adapt(reinterpret_cast<const quaternion_type*>(m_quaternions.data()), {m_quaternions.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::q()
    {
        return xt::adapt(reinterpret_cast<quaternion_type*>(m_quaternions.data()), {m_quaternions.size()});
    }

    // velocity

    template<std::size_t dim>
    auto scopi_container<dim>::v() const
    {
        return xt::adapt(reinterpret_cast<const velocity_type*>(m_velocities.data()), {m_velocities.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::v()
    {
        return xt::adapt(reinterpret_cast<velocity_type*>(m_velocities.data()), {m_velocities.size()});
    }

    // desired velocity

    template<std::size_t dim>
    auto scopi_container<dim>::vd() const
    {
        return xt::adapt(reinterpret_cast<const velocity_type*>(m_desired_velocities.data()), {m_desired_velocities.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::vd()
    {
        return xt::adapt(reinterpret_cast<velocity_type*>(m_desired_velocities.data()), {m_desired_velocities.size()});
    }

    // omega

    template<std::size_t dim>
    auto scopi_container<dim>::omega() const
    {
        return xt::adapt(reinterpret_cast<const rotation_type*>(m_omega.data()), {m_omega.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::omega()
    {
        return xt::adapt(reinterpret_cast<rotation_type*>(m_omega.data()), {m_omega.size()});
    }

    // desired velocity

    template<std::size_t dim>
    auto scopi_container<dim>::desired_omega() const
    {
        return xt::adapt(reinterpret_cast<const rotation_type*>(m_desired_omega.data()), {m_desired_omega.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::desired_omega()
    {
        return xt::adapt(reinterpret_cast<rotation_type*>(m_desired_omega.data()), {m_desired_omega.size()});
    }

    // force

    template<std::size_t dim>
    auto scopi_container<dim>::f() const
    {
        return xt::adapt(reinterpret_cast<const force_type*>(m_forces.data()), {m_forces.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::f()
    {
        return xt::adapt(reinterpret_cast<force_type*>(m_forces.data()), {m_forces.size()});
    }

    // mass

    template<std::size_t dim>
    auto scopi_container<dim>::m() const
    {
        return xt::adapt(reinterpret_cast<const mass_type*>(m_masses.data()), {m_masses.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::m()
    {
        return xt::adapt(reinterpret_cast<mass_type*>(m_masses.data()), {m_masses.size()});
    }

    // moment of inertia

    template<std::size_t dim>
    auto scopi_container<dim>::j() const
    {
        return xt::adapt(reinterpret_cast<const moment_type*>(m_moments_inertia.data()), {m_moments_inertia.size()});
    }

    template<std::size_t dim>
    auto scopi_container<dim>::j()
    {
        return xt::adapt(reinterpret_cast<moment_type*>(m_moments_inertia.data()), {m_moments_inertia.size()});
    }

    template<std::size_t dim>
    void scopi_container<dim>::reset_periodic()
    {
        m_periodic_indices.clear();
        std::size_t size = m_periodic_ptr;
        m_positions.resize(size);
        m_quaternions.resize(size);
        m_velocities.resize(size);
        m_desired_velocities.resize(size);
        m_omega.resize(size);
        m_desired_omega.resize(size);
        m_forces.resize(size);
        m_masses.resize(size);
        m_moments_inertia.resize(size);
        m_offset.resize(size+1);
        m_shapes_id.resize(size);

        m_periodic_added = false;
        m_periodic_ptr = size;
    }

    template<std::size_t dim>
    std::size_t scopi_container<dim>::object_index(std::size_t i) const
    {
        auto lower = std::upper_bound(m_offset.cbegin(), m_offset.cend(), i);
        return std::distance(m_offset.cbegin(), lower)-1;
    }

    template<std::size_t dim>
    std::size_t scopi_container<dim>::offset(std::size_t i) const
    {
        return m_offset[i];
    }

}
