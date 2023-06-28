#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        Vector2D step = (end - start) / (num_nodes - 1);
        for (int i = 0; i < num_nodes; i++) {
            Mass* mass = new Mass(start + step * i, node_mass, false);
            mass->velocity = Vector2D(0, 0);
            if (i > 0) {
                Spring* spring = new Spring(masses.back(), mass, k);
                springs.push_back(spring);
            }
            masses.push_back(mass);
        }

        for (auto& i : pinned_nodes) {
            masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto& s : springs)
        {
            Vector2D dir = (s->m2->position - s->m1->position).unit();
            double dis = (s->m2->position - s->m1->position).norm();
            Vector2D res = s->k * (dis - s->rest_length) * dir;

            s->m1->forces += res;
            s->m2->forces -= res;
        }


        double damping_factor = 8e-4;
        for (auto& m : masses)
        {
            if (!m->pinned)
            {
                m->forces += gravity;
                m->forces -= damping_factor * m->velocity;
                Vector2D a = m->forces / m->mass;
                m->velocity += a * delta_t;
                // m->position += m->velocity * delta_t;
                Vector2D vplus = m->velocity + a * delta_t;
                m->position += vplus * delta_t;
            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
            }
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto& s : springs)
        {
            Vector2D dir = (s->m2->position - s->m1->position).unit();
            double dis = (s->m2->position - s->m1->position).norm();
            Vector2D res = s->k * (dis - s->rest_length) * dir; // s->rest_length is the original length!

            s->m1->forces += res;
            s->m2->forces -= res;
        }

        double damping_factor = 8e-5;
        for (auto& m : masses)
        {
            if (!m->pinned)
            {
                Vector2D a = m->forces / m->mass + gravity;
                // Vector2D next_position = m->position + (m->position - m->last_position) + a * delta_t * delta_t;
                Vector2D next_position = m->position + (1 - damping_factor) * (m->position - m->last_position) + a * delta_t * delta_t;

                m->last_position = m->position;
                m->position = next_position;
            }

            m->forces = Vector2D(0, 0);
        }
    }
}
