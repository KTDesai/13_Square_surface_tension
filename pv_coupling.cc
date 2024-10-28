#include "pv_coupling.hh"
#include<iostream>
#include<fstream>
#include <vector>

PV_coupling::PV_coupling(int no_of_cells, double nu_value, double initial_pressure, double initial_alpha, double initial_rho, double initial_mu, Vector initial_velocity): nu(no_of_cells), velocity(no_of_cells), pressure(no_of_cells), alpha(no_of_cells), rho(no_of_cells), mu(no_of_cells)
{
    nu.set_scalar_field_values(nu_value);
    velocity.set_vector_field_values(initial_velocity);
    pressure.set_scalar_field_values(initial_pressure);
    alpha.set_scalar_field_values(initial_alpha);
    rho.set_scalar_field_values(initial_rho);
    mu.set_scalar_field_values(initial_mu);
    
    velocity_a_matrix_rate_of_change.resize(no_of_cells, no_of_cells);
    velocity_a_matrix_diffusion.resize(no_of_cells, no_of_cells);
    velocity_a_matrix_convection.resize(no_of_cells, no_of_cells);
    velocity_a_matrix_combined.resize(no_of_cells, no_of_cells);
    

    velocity_b_vector_rate_of_change_ux.resize(no_of_cells);
    velocity_b_vector_rate_of_change_uy.resize(no_of_cells);
    velocity_b_vector_source_ux.resize(no_of_cells);
    velocity_b_vector_source_uy.resize(no_of_cells);
    velocity_b_vector_diffusion_ux.resize(no_of_cells);
    velocity_b_vector_diffusion_uy.resize(no_of_cells);
    velocity_b_vector_convection_ux.resize(no_of_cells);
    velocity_b_vector_convection_uy.resize(no_of_cells);
    velocity_b_vector_combined_ux.resize(no_of_cells);
    velocity_b_vector_combined_uy.resize(no_of_cells);

   velocity_a_matrix_rate_of_change.setZero();
   velocity_a_matrix_diffusion.setZero();
   velocity_a_matrix_convection.setZero();
   velocity_a_matrix_combined.setZero();

   velocity_b_vector_rate_of_change_ux.setZero();
   velocity_b_vector_rate_of_change_uy.setZero();
   velocity_b_vector_source_ux.setZero();
   velocity_b_vector_source_uy.setZero();
   velocity_b_vector_diffusion_ux.setZero();
   velocity_b_vector_diffusion_uy.setZero();
   velocity_b_vector_convection_ux.setZero();
   velocity_b_vector_convection_uy.setZero();
   velocity_b_vector_combined_ux.setZero();
   velocity_b_vector_combined_uy.setZero();

    pressure_a_matrix_rate_of_change.resize(no_of_cells, no_of_cells);
    pressure_a_matrix_diffusion.resize(no_of_cells, no_of_cells);
    pressure_a_matrix_convection.resize(no_of_cells, no_of_cells);
    pressure_a_matrix_combined.resize(no_of_cells, no_of_cells);

    pressure_b_vector_rate_of_change.resize(no_of_cells);
    pressure_b_vector_source.resize(no_of_cells);
    pressure_b_vector_diffusion.resize(no_of_cells);
    pressure_b_vector_convection.resize(no_of_cells);
    pressure_b_vector_combined.resize(no_of_cells);

   pressure_a_matrix_rate_of_change.setZero();
   pressure_a_matrix_diffusion.setZero();
   pressure_a_matrix_convection.setZero();
   pressure_a_matrix_combined.setZero();

   pressure_b_vector_rate_of_change.setZero();
   pressure_b_vector_source.setZero();
   pressure_b_vector_diffusion.setZero();
   pressure_b_vector_convection.setZero();
   pressure_b_vector_combined.setZero();

    alpha_a_matrix_rate_of_change.resize(no_of_cells, no_of_cells);
    alpha_a_matrix_convection.resize(no_of_cells, no_of_cells);
    alpha_a_matrix_combined.resize(no_of_cells, no_of_cells);

    alpha_b_vector_rate_of_change.resize(no_of_cells);
    alpha_b_vector_source.resize(no_of_cells);
    alpha_b_vector_convection.resize(no_of_cells);
    alpha_b_vector_combined.resize(no_of_cells);

   alpha_a_matrix_rate_of_change.setZero();
   alpha_a_matrix_convection.setZero();
   alpha_a_matrix_combined.setZero();

   alpha_b_vector_rate_of_change.setZero();
   alpha_b_vector_source.setZero();
   alpha_b_vector_convection.setZero();
   alpha_b_vector_combined.setZero();

}

void PV_coupling::velocity_reset_rate_of_change()
{
    velocity_a_matrix_rate_of_change.setZero();
    velocity_b_vector_rate_of_change_ux.setZero();
    velocity_b_vector_rate_of_change_uy.setZero();
}

void PV_coupling::velocity_reset_source()
{
    velocity_b_vector_source_ux.setZero();
    velocity_b_vector_source_uy.setZero();
}

void PV_coupling::velocity_reset_diffusion()
{
    velocity_a_matrix_diffusion.setZero();
    velocity_b_vector_diffusion_ux.setZero();
    velocity_b_vector_diffusion_uy.setZero();
}

void PV_coupling::velocity_reset_convection()
{
    velocity_a_matrix_convection.setZero();
    velocity_b_vector_convection_ux.setZero();
    velocity_b_vector_convection_uy.setZero();

}

void PV_coupling::velocity_reset_combined()
{
    velocity_a_matrix_combined.setZero();
    velocity_b_vector_combined_ux.setZero();
    velocity_b_vector_combined_uy.setZero();
}

void PV_coupling::velocity_update_a_matrix_rate_of_change(double a, int index)
{
   velocity_a_matrix_rate_of_change(index, index) = velocity_a_matrix_rate_of_change(index, index) + a ;         
}

void PV_coupling::velocity_update_b_vector_rate_of_change_ux(double a, int index)
{
    velocity_b_vector_rate_of_change_ux(index) = velocity_b_vector_rate_of_change_ux(index) + a;
}

void PV_coupling::velocity_update_b_vector_rate_of_change_uy(double a, int index)
{
    velocity_b_vector_rate_of_change_uy(index) = velocity_b_vector_rate_of_change_uy(index) + a;
}

void PV_coupling::velocity_update_b_vector_source_ux(double a, int index)
{
    velocity_b_vector_source_ux(index) = velocity_b_vector_source_ux(index) + a;
}

void PV_coupling::velocity_update_b_vector_source_uy(double a, int index)
{
    velocity_b_vector_source_uy(index) = velocity_b_vector_source_uy(index) + a;
}

void PV_coupling::velocity_update_a_matrix_diffusion(double a, int row_index, int column_index)
{ 
    velocity_a_matrix_diffusion(row_index, column_index) = velocity_a_matrix_diffusion(row_index, column_index) + a;                
}                                                                         

void PV_coupling::velocity_update_b_vector_diffusion_ux(double a, int index)
{
    velocity_b_vector_diffusion_ux(index) = velocity_b_vector_diffusion_ux(index) + a;           
}

void PV_coupling::velocity_update_b_vector_diffusion_uy(double a, int index)
{
    velocity_b_vector_diffusion_uy(index) = velocity_b_vector_diffusion_uy(index) + a;           
}

void PV_coupling::velocity_update_a_matrix_convection(double a, int row_index, int column_index)
{ 
    velocity_a_matrix_convection(row_index, column_index) = velocity_a_matrix_convection(row_index, column_index) + a;                 
}                                                                         

void PV_coupling::velocity_update_b_vector_convection_ux(double a, int index)
{
    velocity_b_vector_convection_ux(index) = velocity_b_vector_convection_ux(index) + a;           
}

void PV_coupling::velocity_update_b_vector_convection_uy(double a, int index)
{
    velocity_b_vector_convection_uy(index) = velocity_b_vector_convection_uy(index) + a;           
}

void PV_coupling::velocity_combine_a_matrices()
{
    for(int i = 0; i<velocity_a_matrix_combined.rows(); i++)
    {
        for(int j=0; j<velocity_a_matrix_combined.cols(); j++)
        {
            velocity_a_matrix_combined(i,j) = velocity_a_matrix_rate_of_change(i,j) + velocity_a_matrix_diffusion(i,j) + velocity_a_matrix_convection(i,j);
        }
    }

}

void PV_coupling::velocity_combine_b_matrices()
{
    for(int k =0 ; k <velocity_b_vector_combined_ux.size(); k++)
    {
        velocity_b_vector_combined_ux(k) = velocity_b_vector_rate_of_change_ux(k) + velocity_b_vector_source_ux(k) + velocity_b_vector_diffusion_ux(k) + velocity_b_vector_convection_ux(k) ;
    }

    for(int k =0 ; k <velocity_b_vector_combined_ux.size(); k++)
    {
        velocity_b_vector_combined_uy(k) = velocity_b_vector_rate_of_change_uy(k) + velocity_b_vector_source_uy(k) + velocity_b_vector_diffusion_uy(k) + velocity_b_vector_convection_uy(k) ;
    }
}

 void PV_coupling::velocity_compute_rate_of_change_matrix(std::vector<Cell> list_of_cells, double delta_time)                                                               
 {
      velocity_reset_rate_of_change();  

      for(int i = 0; i<list_of_cells.size(); i++)
      {  
            int cell_index = list_of_cells[i].get_cell_index();
            double cell_vol = list_of_cells[cell_index].get_cell_volume();
            double velocity_x_old = velocity.vector_field_values_x(cell_index);
            double velocity_y_old = velocity.vector_field_values_y(cell_index);
            double rho_new = rho.scalar_field_values(cell_index);
            double rho_old = rho.scalar_old_field_values(cell_index);

            double rate_of_change_diagonal_contribution = (cell_vol * rho_new)/delta_time;
            double rate_of_change_source_contribution_ux = (cell_vol * rho_old * velocity_x_old)/ delta_time;
            double rate_of_change_source_contribution_uy = (cell_vol * rho_old * velocity_y_old)/ delta_time;

             velocity_update_a_matrix_rate_of_change(rate_of_change_diagonal_contribution, cell_index);
             velocity_update_b_vector_rate_of_change_ux(rate_of_change_source_contribution_ux, cell_index);
             velocity_update_b_vector_rate_of_change_uy(rate_of_change_source_contribution_uy, cell_index);
      }
 } 


void PV_coupling::velocity_under_relaxation(std::vector<Cell> list_of_cells, double alpha)                                                               
 {  
    double alpha_source_factor = (1 - alpha)/alpha ;

    for(int i=0; i<list_of_cells.size(); i++)
    {
         velocity_b_vector_combined_ux(i)  =  velocity_b_vector_combined_ux(i) + (alpha_source_factor* velocity_a_matrix_combined(i,i) * velocity.vector_field_values_x(i)) ;
         velocity_b_vector_combined_uy(i)  =  velocity_b_vector_combined_uy(i) + (alpha_source_factor* velocity_a_matrix_combined(i,i) * velocity.vector_field_values_y(i)) ; 
         velocity_a_matrix_combined(i,i) =  velocity_a_matrix_combined(i,i) / alpha;
    }
 }  


void PV_coupling::velocity_compute_diffusion_matrix_old(std::vector<Face> list_of_faces, std::vector<Boundary> list_of_boundaries, Vector_boundary_field velocity_boundary)                                                               
 { 
      velocity_reset_diffusion(); 

       for(int i = 0; i<list_of_faces.size(); i++)
        {  
          double a_owner_diagonal, a_neighbour_diagonal, a_owner_neighbour, a_neighbour_owner, b_owner_source_ux, b_owner_source_uy, b_neighbour_source_ux, b_neighbour_source_uy;

          int face_owner_index = list_of_faces[i].get_owner_index();
          int face_neighbour_index = list_of_faces[i].get_neighbour_index();

          double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
          double face_delta = list_of_faces[i].get_face_delta(); 

          Vector face_normal = list_of_faces[i].get_face_normal();

          double face_normal_magnitude = face_normal.magnitude_of_vector();

          if(list_of_faces[i].get_is_boundary_face() == true)
            {       
              int boundary_index = list_of_faces[i].get_boundary_index();
              Boundary b = list_of_boundaries[boundary_index];

                //if(b.return_boundary_name() == "empty")
                if(velocity_boundary.get_vector_boundary_name(boundary_index) == "empty")
                  {
                    continue;
                  }

                else
                  {
                    double nu_owner = nu.scalar_field_values(face_owner_index);
                    double nu_boundary_face = interpolation_factor * nu_owner;

                    //if(b.get_is_vector_dirichlet() == true)
                    if(velocity_boundary.get_vector_boundary_type(boundary_index) == "dirichlet")
                    {   
                        // std::vector<double> dirichlet_value = b.get_vector_dirichlet_boundary_value();
                        std::vector<std::vector<double>> dirichlet_value = velocity_boundary.get_vector_boundary_field_values();


                        a_owner_diagonal = nu_boundary_face*face_normal_magnitude*face_delta;
                        velocity_update_a_matrix_diffusion(a_owner_diagonal, face_owner_index, face_owner_index); 

                        b_owner_source_ux = nu_boundary_face*face_normal_magnitude*face_delta*dirichlet_value[boundary_index][0];
                        velocity_update_b_vector_diffusion_ux(b_owner_source_ux, face_owner_index) ;

                        b_owner_source_uy = nu_boundary_face*face_normal_magnitude*face_delta*dirichlet_value[boundary_index][1];
                        velocity_update_b_vector_diffusion_uy(b_owner_source_uy, face_owner_index) ;

                    }

                    else
                    {
                        // std::vector<double> neumann_gradient =  b.get_vector_neumann_boundary_value();
                        std::vector<std::vector<double>> neumann_gradient = velocity_boundary.get_vector_boundary_field_values();

                        b_owner_source_ux = nu_owner * face_normal_magnitude * neumann_gradient[boundary_index][0];
                        velocity_update_b_vector_diffusion_ux(b_owner_source_ux, face_owner_index);

                        b_owner_source_uy = nu_owner * face_normal_magnitude * neumann_gradient[boundary_index][1];
                        velocity_update_b_vector_diffusion_uy(b_owner_source_ux, face_owner_index);
                    }
                  }
            }

                else
                {
                   double nu_owner_cell = nu.scalar_field_values(face_owner_index);
                   double nu_neighbour_cell = nu.scalar_field_values(face_neighbour_index);

                   double nu_face = (interpolation_factor*nu_owner_cell) + ((1 - interpolation_factor)*(nu_neighbour_cell)) ;

                   a_owner_diagonal = nu_face * face_normal_magnitude * face_delta;
                   a_neighbour_diagonal = nu_face * face_normal_magnitude * face_delta;

                   a_owner_neighbour =  -1 * nu_face * face_normal_magnitude * face_delta;
                   a_neighbour_owner =  -1 * nu_face * face_normal_magnitude * face_delta;
                   
                   velocity_update_a_matrix_diffusion(a_owner_diagonal, face_owner_index, face_owner_index);
                   velocity_update_a_matrix_diffusion(a_neighbour_diagonal, face_neighbour_index, face_neighbour_index);
                   velocity_update_a_matrix_diffusion(a_owner_neighbour, face_owner_index, face_neighbour_index);
                   velocity_update_a_matrix_diffusion(a_neighbour_owner, face_neighbour_index, face_owner_index);
            
                }
        }     
 }

void PV_coupling::velocity_compute_diffusion_matrix(std::vector<Face> list_of_faces, std::vector<Boundary> list_of_boundaries, Vector_boundary_field velocity_boundary)                                                               
 { 
      velocity_reset_diffusion(); 

       for(int i = 0; i<list_of_faces.size(); i++)
        {  
          double a_owner_diagonal, a_neighbour_diagonal, a_owner_neighbour, a_neighbour_owner, b_owner_source_ux, b_owner_source_uy, b_neighbour_source_ux, b_neighbour_source_uy;

          int face_owner_index = list_of_faces[i].get_owner_index();
          int face_neighbour_index = list_of_faces[i].get_neighbour_index();

          double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
          double face_delta = list_of_faces[i].get_face_delta(); 

          Vector face_normal = list_of_faces[i].get_face_normal();

          double face_normal_magnitude = face_normal.magnitude_of_vector();

          if(list_of_faces[i].get_is_boundary_face() == true)
            {       
              int boundary_index = list_of_faces[i].get_boundary_index();
              Boundary b = list_of_boundaries[boundary_index];

                //if(b.return_boundary_name() == "empty")
                if(velocity_boundary.get_vector_boundary_name(boundary_index) == "empty")
                  {
                    continue;
                  }

                else
                  {
                    double mu_owner = mu.scalar_field_values(face_owner_index);
                    double mu_boundary_face = interpolation_factor * mu_owner;

                    //if(b.get_is_vector_dirichlet() == true)
                    if(velocity_boundary.get_vector_boundary_type(boundary_index) == "dirichlet")
                    {   
                        // std::vector<double> dirichlet_value = b.get_vector_dirichlet_boundary_value();
                        std::vector<std::vector<double>> dirichlet_value = velocity_boundary.get_vector_boundary_field_values();


                        a_owner_diagonal = mu_boundary_face*face_normal_magnitude*face_delta;
                        velocity_update_a_matrix_diffusion(a_owner_diagonal, face_owner_index, face_owner_index); 

                        b_owner_source_ux = mu_boundary_face*face_normal_magnitude*face_delta*dirichlet_value[boundary_index][0];
                        velocity_update_b_vector_diffusion_ux(b_owner_source_ux, face_owner_index) ;

                        b_owner_source_uy = mu_boundary_face*face_normal_magnitude*face_delta*dirichlet_value[boundary_index][1];
                        velocity_update_b_vector_diffusion_uy(b_owner_source_uy, face_owner_index) ;

                    }

                    else
                    {
                        // std::vector<double> neumann_gradient =  b.get_vector_neumann_boundary_value();
                        std::vector<std::vector<double>> neumann_gradient = velocity_boundary.get_vector_boundary_field_values();

                        b_owner_source_ux = mu_owner * face_normal_magnitude * neumann_gradient[boundary_index][0];
                        velocity_update_b_vector_diffusion_ux(b_owner_source_ux, face_owner_index);

                        b_owner_source_uy = mu_owner * face_normal_magnitude * neumann_gradient[boundary_index][1];
                        velocity_update_b_vector_diffusion_uy(b_owner_source_ux, face_owner_index);
                    }
                  }
            }

                else
                {
                   double mu_owner_cell = mu.scalar_field_values(face_owner_index);
                   double mu_neighbour_cell = mu.scalar_field_values(face_neighbour_index);

                   double mu_face = (interpolation_factor*mu_owner_cell) + ((1 - interpolation_factor)*(mu_neighbour_cell)) ;

                   a_owner_diagonal = mu_face * face_normal_magnitude * face_delta;
                   a_neighbour_diagonal = mu_face * face_normal_magnitude * face_delta;

                   a_owner_neighbour =  -1 * mu_face * face_normal_magnitude * face_delta;
                   a_neighbour_owner =  -1 * mu_face * face_normal_magnitude * face_delta;
                   
                   velocity_update_a_matrix_diffusion(a_owner_diagonal, face_owner_index, face_owner_index);
                   velocity_update_a_matrix_diffusion(a_neighbour_diagonal, face_neighbour_index, face_neighbour_index);
                   velocity_update_a_matrix_diffusion(a_owner_neighbour, face_owner_index, face_neighbour_index);
                   velocity_update_a_matrix_diffusion(a_neighbour_owner, face_neighbour_index, face_owner_index);
            
                }
        }     
 }


 void PV_coupling::velocity_compute_convection_matrix_old(std::vector<Face> list_of_faces, std::vector<Boundary> list_of_boundaries, Vector_boundary_field velocity_boundary)                                                               
 { 
      velocity_reset_convection();  

      for(int i = 0; i<list_of_faces.size(); i++)
        {   
             Vector face_velocity;
             int face_index = list_of_faces[i].get_face_index();

             double a_owner_diagonal, a_neighbour_diagonal, a_owner_neighbour, a_neighbour_owner, b_owner_source_ux, b_neighbour_source_ux, b_owner_source_uy, b_neighbour_source_uy;
             double neighbour_velocity_x, neighbour_velocity_y, neighbour_velocity_z;

             int owner_index = list_of_faces[i].get_owner_index();
             int neighbour_index = list_of_faces[i].get_neighbour_index();

             double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
             double face_delta = list_of_faces[i].get_face_interpolation_factor();

            if(list_of_faces[i].get_is_boundary_face() == true)
            {
                int boundary_index = list_of_faces[i].get_boundary_index();
                Boundary b = list_of_boundaries[list_of_faces[i].get_boundary_index()];
                // std::vector<double> vel = b.get_vector_dirichlet_boundary_value();
                std::vector<std::vector<double>> vel = velocity_boundary.get_vector_boundary_field_values();

                face_velocity = Vector(vel[boundary_index][0], vel[boundary_index][1], vel[boundary_index][2]);
            }

            else
            {      
                double owner_vx = velocity.vector_field_values_x(owner_index);
                double owner_vy = velocity.vector_field_values_y(owner_index);
                double owner_vz = velocity.vector_field_values_z(owner_index);
                double neighbour_vx = velocity.vector_field_values_x(neighbour_index);
                double neighbour_vy = velocity.vector_field_values_y(neighbour_index);
                double neighbour_vz = velocity.vector_field_values_z(neighbour_index);
                Vector owner_velocity(owner_vx, owner_vy, owner_vz);
                Vector neighbour_velocity(neighbour_vx, neighbour_vy, neighbour_vz);
                face_velocity = owner_velocity * interpolation_factor + neighbour_velocity * (1 - interpolation_factor);
            }

                Vector face_normal = list_of_faces[i].get_face_normal();

                double velocity_flux = face_normal.dot_product(face_velocity);

                double face_normal_magnitude = face_normal.magnitude_of_vector();

                if(velocity_flux >= 0)
                {
                    interpolation_factor = 1;
                }

                else
                {
                    interpolation_factor = 0;
                }

                 if(list_of_faces[i].get_is_boundary_face() == true)
                {
                     
                     int boundary_index = list_of_faces[i].get_boundary_index();
                     Boundary b = list_of_boundaries[boundary_index];

                    //  if(b.return_boundary_name() == "empty")
                    if(velocity_boundary.get_vector_boundary_name(boundary_index) == "empty")
                     {
                        continue;
                     }

                     else
                     {
                        //if(b.get_is_vector_dirichlet() == true)
                        if(velocity_boundary.get_vector_boundary_type(boundary_index) == "dirichlet")
                        {   
                            // std::vector<double> dirichlet_value = b.get_vector_dirichlet_boundary_value();
                            std::vector<std::vector<double>> dirichlet_value = velocity_boundary.get_vector_boundary_field_values();

                            b_owner_source_ux = -1*velocity_flux * dirichlet_value[boundary_index][0];
                            velocity_update_b_vector_convection_ux(b_owner_source_ux, owner_index);

                            b_owner_source_uy = -1*velocity_flux * dirichlet_value[boundary_index][1];
                            velocity_update_b_vector_convection_uy(b_owner_source_uy, owner_index);
                        }

                       //if(b.get_is_vector_neumann() == true)
                       if(velocity_boundary.get_vector_boundary_type(boundary_index) == "neumann")
                        {
                            //std::vector<double> neumann_gradient =  b.get_vector_neumann_boundary_value();
                            std::vector<std::vector<double>> neumann_gradient = velocity_boundary.get_vector_boundary_field_values();

                            a_owner_diagonal = velocity_flux;
                            velocity_update_a_matrix_convection(a_owner_diagonal, owner_index, owner_index);

                            b_owner_source_ux = -1*(velocity_flux/face_delta) * neumann_gradient[boundary_index][0];
                            velocity_update_b_vector_convection_ux(b_owner_source_ux, owner_index);

                            b_owner_source_uy = -1*(velocity_flux/face_delta) * neumann_gradient[boundary_index][1];
                            velocity_update_b_vector_convection_uy(b_owner_source_uy, owner_index);
                        }
                     }
                }

                else
                {       

                  a_owner_diagonal = interpolation_factor*velocity_flux;
                  a_neighbour_diagonal = -1*(1 - interpolation_factor)*velocity_flux;

                   velocity_update_a_matrix_convection(a_owner_diagonal, owner_index, owner_index);
                   velocity_update_a_matrix_convection(a_neighbour_diagonal, neighbour_index, neighbour_index);

                   a_owner_neighbour = (1 - interpolation_factor)*velocity_flux;
                   a_neighbour_owner = -1*velocity_flux*interpolation_factor;
                   
            
                   velocity_update_a_matrix_convection(a_owner_neighbour, owner_index, neighbour_index);
                   velocity_update_a_matrix_convection(a_neighbour_owner, neighbour_index, owner_index);
                   
                } 
        }
    }  

void PV_coupling::velocity_compute_convection_matrix(std::vector<Face> list_of_faces, std::vector<Boundary> list_of_boundaries, Vector_boundary_field velocity_boundary)                                                               
 { 
      velocity_reset_convection();  

      for(int i = 0; i<list_of_faces.size(); i++)
        {   
             Vector face_velocity;
             int face_index = list_of_faces[i].get_face_index();

             double a_owner_diagonal, a_neighbour_diagonal, a_owner_neighbour, a_neighbour_owner, b_owner_source_ux, b_neighbour_source_ux, b_owner_source_uy, b_neighbour_source_uy;
             double neighbour_velocity_x, neighbour_velocity_y, neighbour_velocity_z;
             double rho_owner, rho_neighbour, rho_face;

             int owner_index = list_of_faces[i].get_owner_index();
             int neighbour_index = list_of_faces[i].get_neighbour_index();

             double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
             double face_delta = list_of_faces[i].get_face_interpolation_factor();

            if(list_of_faces[i].get_is_boundary_face() == true)
            {
                int boundary_index = list_of_faces[i].get_boundary_index();
                Boundary b = list_of_boundaries[list_of_faces[i].get_boundary_index()];
                // std::vector<double> vel = b.get_vector_dirichlet_boundary_value();
                std::vector<std::vector<double>> vel = velocity_boundary.get_vector_boundary_field_values();

                face_velocity = Vector(vel[boundary_index][0], vel[boundary_index][1], vel[boundary_index][2]);

                rho_owner = rho.scalar_old_field_values (owner_index);
                rho_face = rho_owner * interpolation_factor; 
            }

            else
            {      
                double owner_vx = velocity.vector_field_values_x(owner_index);
                double owner_vy = velocity.vector_field_values_y(owner_index);
                double owner_vz = velocity.vector_field_values_z(owner_index);
                double neighbour_vx = velocity.vector_field_values_x(neighbour_index);
                double neighbour_vy = velocity.vector_field_values_y(neighbour_index);
                double neighbour_vz = velocity.vector_field_values_z(neighbour_index);
                Vector owner_velocity(owner_vx, owner_vy, owner_vz);
                Vector neighbour_velocity(neighbour_vx, neighbour_vy, neighbour_vz);
                face_velocity = owner_velocity * interpolation_factor + neighbour_velocity * (1 - interpolation_factor);

                rho_owner = rho.scalar_old_field_values (owner_index);
                rho_neighbour = rho.scalar_old_field_values(neighbour_index);

                rho_face = (interpolation_factor * rho_owner) + ((1.0 - interpolation_factor) * rho_neighbour);
            }

                Vector face_normal = list_of_faces[i].get_face_normal();

                double velocity_flux = face_normal.dot_product(face_velocity);

                velocity_flux = velocity_flux * rho_face;

                double face_normal_magnitude = face_normal.magnitude_of_vector();

                if(velocity_flux >= 0)
                {
                    interpolation_factor = 1;
                }

                else
                {
                    interpolation_factor = 0;
                }

                 if(list_of_faces[i].get_is_boundary_face() == true)
                {
                     
                     int boundary_index = list_of_faces[i].get_boundary_index();
                     Boundary b = list_of_boundaries[boundary_index];

                    //  if(b.return_boundary_name() == "empty")
                    if(velocity_boundary.get_vector_boundary_name(boundary_index) == "empty")
                     {
                        continue;
                     }

                     else
                     {
                        //if(b.get_is_vector_dirichlet() == true)
                        if(velocity_boundary.get_vector_boundary_type(boundary_index) == "dirichlet")
                        {   
                            // std::vector<double> dirichlet_value = b.get_vector_dirichlet_boundary_value();
                            std::vector<std::vector<double>> dirichlet_value = velocity_boundary.get_vector_boundary_field_values();

                            b_owner_source_ux = -1*velocity_flux * dirichlet_value[boundary_index][0];
                            velocity_update_b_vector_convection_ux(b_owner_source_ux, owner_index);

                            b_owner_source_uy = -1*velocity_flux * dirichlet_value[boundary_index][1];
                            velocity_update_b_vector_convection_uy(b_owner_source_uy, owner_index);
                        }

                       //if(b.get_is_vector_neumann() == true)
                       if(velocity_boundary.get_vector_boundary_type(boundary_index) == "neumann")
                        {
                            //std::vector<double> neumann_gradient =  b.get_vector_neumann_boundary_value();
                            std::vector<std::vector<double>> neumann_gradient = velocity_boundary.get_vector_boundary_field_values();

                            a_owner_diagonal = velocity_flux;
                            velocity_update_a_matrix_convection(a_owner_diagonal, owner_index, owner_index);

                            b_owner_source_ux = -1*(velocity_flux/face_delta) * neumann_gradient[boundary_index][0];
                            velocity_update_b_vector_convection_ux(b_owner_source_ux, owner_index);

                            b_owner_source_uy = -1*(velocity_flux/face_delta) * neumann_gradient[boundary_index][1];
                            velocity_update_b_vector_convection_uy(b_owner_source_uy, owner_index);
                        }
                     }
                }

                else
                {       

                  a_owner_diagonal = interpolation_factor*velocity_flux;
                  a_neighbour_diagonal = -1*(1 - interpolation_factor)*velocity_flux;

                   velocity_update_a_matrix_convection(a_owner_diagonal, owner_index, owner_index);
                   velocity_update_a_matrix_convection(a_neighbour_diagonal, neighbour_index, neighbour_index);

                   a_owner_neighbour = (1 - interpolation_factor)*velocity_flux;
                   a_neighbour_owner = -1*velocity_flux*interpolation_factor;
                   
            
                   velocity_update_a_matrix_convection(a_owner_neighbour, owner_index, neighbour_index);
                   velocity_update_a_matrix_convection(a_neighbour_owner, neighbour_index, owner_index);
                   
                } 
        }
    }  


void PV_coupling::velocity_calculate_initial_residuals(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_initial_residual_ux, 
                                               std::ofstream &sum_initial_residual_uy)
{ 
   if(n == 199)
  { 
   std::ofstream file_initial_residuals_ux("initial_residuals_ux.txt");
   std::ofstream file_initial_residuals_uy("initial_residuals_uy.txt");

   velocity_initial_residuals_ux = velocity_b_vector_combined_ux - velocity_a_matrix_combined * velocity.vector_field_values_x;
   velocity_initial_residuals_uy = velocity_b_vector_combined_uy - velocity_a_matrix_combined * velocity.vector_field_values_y;

   for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
            file_initial_residuals_ux<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<velocity_initial_residuals_ux(j + (i*x_distance.size()))<<std::endl;
            file_initial_residuals_uy<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<velocity_initial_residuals_uy(j + (i*x_distance.size()))<<std::endl;
      }
   } 

   file_initial_residuals_ux.close();
   file_initial_residuals_uy.close();
  } 
 
  velocity_initial_residuals_ux = velocity_b_vector_combined_ux - velocity_a_matrix_combined * velocity.vector_field_values_x;
  velocity_initial_residuals_uy = velocity_b_vector_combined_uy - velocity_a_matrix_combined * velocity.vector_field_values_y;

  double sum_r_ux = 0.0, sum_r_uy = 0.0;

  for(int i=0; i<velocity_initial_residuals_ux.size(); i++)
  {
     sum_r_ux = sum_r_ux + velocity_initial_residuals_ux(i);
     sum_r_uy = sum_r_uy + velocity_initial_residuals_uy(i);
  } 

  sum_initial_residual_ux<<n<<"\t"<<sum_r_ux<<std::endl;
  sum_initial_residual_uy<<n<<"\t"<<sum_r_ux<<std::endl;
 
}    


void PV_coupling::velocity_solve_matrices()
{
    // velocity.vector_field_values_x = velocity_a_matrix_combined.lu().solve(velocity_b_vector_combined_ux);
    // velocity.vector_field_values_y = velocity_a_matrix_combined.lu().solve(velocity_b_vector_combined_uy);

    velocity.vector_field_values_x = velocity_a_matrix_combined.fullPivHouseholderQr().solve(velocity_b_vector_combined_ux);
    velocity.vector_field_values_y = velocity_a_matrix_combined.fullPivHouseholderQr().solve(velocity_b_vector_combined_uy);

    velocity.set_old_vector_field_values();
}  


void PV_coupling::velocity_calculate_final_residuals(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_final_residual_ux, std::ofstream &sum_final_residual_uy)
{
   if(n==199)
 {  
   std::ofstream file_final_residuals_ux("final_residuals_ux.txt");
   std::ofstream file_final_residuals_uy("final_residuals_uy.txt");
   
   velocity_final_residuals_ux = velocity_b_vector_combined_ux - velocity_a_matrix_combined * velocity.vector_field_values_x;
   velocity_final_residuals_uy = velocity_b_vector_combined_uy - velocity_a_matrix_combined * velocity.vector_field_values_y;

      for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
        file_final_residuals_ux<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<velocity_final_residuals_ux(j + (i*x_distance.size()))<<std::endl;
        file_final_residuals_uy<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<velocity_final_residuals_uy(j + (i*x_distance.size()))<<std::endl;
      }

   } 

   file_final_residuals_ux.close();
   file_final_residuals_uy.close();
 }

  velocity_final_residuals_ux = velocity_b_vector_combined_ux - velocity_a_matrix_combined * velocity.vector_field_values_x;
  velocity_final_residuals_uy = velocity_b_vector_combined_uy - velocity_a_matrix_combined * velocity.vector_field_values_y;

  double sum_r_ux = 0.0, sum_r_uy = 0.0;

  for(int i=0; i<velocity_initial_residuals_ux.size(); i++)
  {
     sum_r_ux = sum_r_ux + velocity_final_residuals_ux(i);
     sum_r_uy = sum_r_uy + velocity_final_residuals_uy(i);
  } 

  sum_final_residual_ux<<n<<"\t"<<sum_r_ux<<std::endl;
  sum_final_residual_uy<<n<<"\t"<<sum_r_ux<<std::endl;
}   

void PV_coupling::velocity_set_face_and_cell_fluxes(std::vector<Cell> &list_of_cells, std::vector<Face> &list_of_faces, std::vector<Boundary> list_of_boundaries, Vector_boundary_field velocity_boundary)
{
    std::vector<double> cell_fluxes;

        for(int i=0; i<list_of_faces.size(); i++)
        { 
            double face_flux;
            Vector face_velocity;

            int owner_index = list_of_faces[i].get_owner_index();
            int neighbour_index = list_of_faces[i].get_neighbour_index();
            double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
            
            Vector face_normal = list_of_faces[i].get_face_normal();

            if(list_of_faces[i].get_is_boundary_face() == true)
            { 
                int boundary_index = list_of_faces[i].get_boundary_index();
                Boundary b = list_of_boundaries[boundary_index];
                
                //if(b.return_boundary_name() == "empty")
                if(velocity_boundary.get_vector_boundary_name(boundary_index) == "empty")
                {
                    break;
                }

                //std::vector<double> bound_values = b.get_vector_dirichlet_boundary_value();
                std::vector<std::vector<double>> dirichlet_value = velocity_boundary.get_vector_boundary_field_values();
                face_velocity = Vector(dirichlet_value[boundary_index][0], dirichlet_value[boundary_index][1], dirichlet_value[boundary_index][2]);  
                //face_velocity = Vector(bound_values[0], bound_values[1], bound_values[2]); 
            }  

            else
            { 
                double x_velocity_owner = velocity.vector_field_values_x(owner_index);
                double y_velocity_owner = velocity.vector_field_values_y(owner_index);
                 double z_velocity_owner = velocity.vector_field_values_z(owner_index);
                //double z_velocity_owner = 0.0;

                double x_velocity_neighbour = velocity.vector_field_values_x(neighbour_index);
                double y_velocity_neighbour = velocity.vector_field_values_y(neighbour_index);
                 double z_velocity_neighbour = velocity.vector_field_values_z(neighbour_index);
                //double z_velocity_neighbour = 0.0;

                Vector owner_vel(x_velocity_owner, y_velocity_owner, z_velocity_owner);
                Vector neighbour_vel(x_velocity_neighbour, y_velocity_neighbour, z_velocity_neighbour);

                face_velocity = (owner_vel * interpolation_factor) + (neighbour_vel * (1 - interpolation_factor));
            }  

            face_flux = face_normal.dot_product(face_velocity);  

            list_of_faces[i].set_face_flux(face_flux);

        }
}


void PV_coupling::velocity_compute_source_matrix(std::vector<Cell> list_of_cells, std::vector<Face> list_of_faces)                                                               
 {  
    velocity_reset_source();

    Vector_field alpha_gauss_grad = alpha.compute_gauss_gradient(list_of_cells, list_of_faces); 
    double sigma  = 0.07;

    // for(int i = 0; i<list_of_cells.size(); i++)
    // {  
    //     int cell_index = list_of_cells[i].get_cell_index();

    //     double cell_volume = list_of_cells[i].get_cell_volume();
    //     double gravity_constant = -9.81; 
    //     double cell_density = rho.scalar_field_values(cell_index);

    //     double gravity_source_term_ux = 0.0;
    //     double gravity_source_term_uy = cell_density * gravity_constant * cell_volume;

    //     velocity_update_b_vector_source_ux(gravity_source_term_ux, cell_index);
    //     velocity_update_b_vector_source_uy(gravity_source_term_uy, cell_index);
    // }    

     for(int i=0; i<list_of_faces.size(); i++)
     {
       
       if(list_of_faces[i].get_is_boundary_face() == true)
       {
          continue;
       }
       
       else
       {
         int owner_index = list_of_faces[i].get_owner_index();
         int neighbour_index = list_of_faces[i].get_neighbour_index();

         Vector face_normal = list_of_faces[i].get_face_normal();

         Vector owner_cell_centre = list_of_cells[owner_index].get_cell_centre();
         Vector neighbour_cell_centre = list_of_cells[neighbour_index].get_cell_centre();

       Vector cell_to_cell_distance = neighbour_cell_centre - owner_cell_centre;
       double magnitude_of_vector_cell_to_cell = cell_to_cell_distance.magnitude_of_vector(); 
  
        Vector cell_to_cell_unit_vector = cell_to_cell_distance / magnitude_of_vector_cell_to_cell;

       

       double alpha_owner_cell = alpha.scalar_field_values(owner_index);
       double alpha_neighbour_cell = alpha.scalar_field_values(neighbour_index);

      double alpha_diff = (alpha_neighbour_cell - alpha_owner_cell)/magnitude_of_vector_cell_to_cell;

      Vector alpha_diff_vector = cell_to_cell_unit_vector * alpha_diff;
      double mag_alpha_diff_vector = alpha_diff_vector.magnitude_of_vector();

      Vector V(0,0,0);
      Vector kappa;

      if(mag_alpha_diff_vector == 0)
      {
         kappa = V;
      }

      else
      {
         kappa = alpha_diff_vector / mag_alpha_diff_vector ;
      }
      
       
      double kappa_face = kappa.dot_product(face_normal);

      double  kappa_face_times_sigma = kappa_face * sigma;

      Vector grad_alpha_owner_cell = Vector(alpha_gauss_grad.vector_field_values_x(owner_index), alpha_gauss_grad.vector_field_values_y(owner_index), alpha_gauss_grad.vector_field_values_z(owner_index));
      Vector grad_alpha_neighbour_cell = Vector(alpha_gauss_grad.vector_field_values_x(neighbour_index), alpha_gauss_grad.vector_field_values_y(neighbour_index), alpha_gauss_grad.vector_field_values_z(neighbour_index));

      Vector f_sigma_owner = grad_alpha_owner_cell * kappa_face_times_sigma;
      Vector f_sigma_neighbour = grad_alpha_neighbour_cell * kappa_face_times_sigma;

      double f_sigma_ux_owner = f_sigma_owner[0];
      double f_sigma_uy_owner = f_sigma_owner[1];
      double f_sigma_uz_owner = f_sigma_owner[2];

      velocity_update_b_vector_source_ux(f_sigma_ux_owner, owner_index);
      velocity_update_b_vector_source_uy(f_sigma_uy_owner, owner_index);

      double f_sigma_ux_neighbour = f_sigma_neighbour[0];
      f_sigma_ux_neighbour = f_sigma_ux_neighbour * -1;

      double f_sigma_uy_neighbour = f_sigma_neighbour[1];
      f_sigma_uy_neighbour = f_sigma_uy_neighbour * -1;

      double f_sigma_uz_neighbour = f_sigma_neighbour[2];
      f_sigma_uz_neighbour = f_sigma_uz_neighbour * -1;

      velocity_update_b_vector_source_ux(f_sigma_ux_neighbour, neighbour_index);
      velocity_update_b_vector_source_uy(f_sigma_uy_neighbour, neighbour_index);

      }
       
     }
 } 

std::vector<double> PV_coupling::velocity_store_ap_coefficients()
{
    std::vector<double> diag_coeffs(velocity_b_vector_combined_ux.size());

    for(int i=0; i<velocity_b_vector_combined_ux.size(); i++)
    {
        diag_coeffs[i] = 1.0/velocity_a_matrix_combined(i,i);
    }

    return diag_coeffs;
}


  void PV_coupling::velocity_correct_cell_centre_velocities(std::vector<Cell> &list_of_cells, Vector_field grad_p , std::vector<double> temp_ap_coeffs)
  {
     for(int i=0; i<list_of_cells.size(); i++)
    {
         velocity.vector_field_values_x(i) -= temp_ap_coeffs[i] * grad_p.vector_field_values_x(i);
         velocity.vector_field_values_y(i) -= temp_ap_coeffs[i] * grad_p.vector_field_values_y(i);
    }
  }

  void PV_coupling::velocity_compute_div_u(std::vector<Cell> list_of_cells, int n_cells, std::vector<Face> list_of_faces)                                                               
 {  
    std::vector<double> div_u;
    div_u.resize(n_cells);
    
    for(int j=0; j<list_of_cells.size(); j++)
    {
      div_u[j] = 0.0;
    }

            //   std::cout<<"In the in the in the in the int loop"<<std::endl;

            // for(int m =0 ; m<list_of_faces.size(); m++)
            // {
            //     std::cout<<"Innnnnnn the loop"<<list_of_faces[m].get_face_flux()<<std::endl;
            // }
            

//  for(int i = 0; i<list_of_cells.size(); i++)
//       {           
  
//             double sum_face_flux = 0.0, face_flux;

//             int cell_index = list_of_cells[i].get_cell_index();

//             std::vector<Face*> list_of_faces  =  list_of_cells[cell_index].get_list_of_faces(); 
           
//             for(int k = 0; k<list_of_faces.size() ; k++)
//             {
//                 Face *f =  list_of_faces[k];

//                 int owner_index = f->get_owner_index();
//                 int neighbour_index = f->get_neighbour_index();

//                 double interpolation_factor = f->get_face_interpolation_factor();
//                 double face_delta = f->get_face_delta();
//                 double face_area = f->get_face_area();

//                 int face_index = f->get_face_index();
               
//                 Vector face_normal = f->get_face_normal();

//                 if(cell_index == neighbour_index)
//                 {
//                     face_normal = face_normal * (-1);
//                 }

//                 double face_normal_magnitude = face_normal.magnitude_of_vector();

//                //face_flux = f->get_face_flux();
//                face_flux = list_of_faces[face_index].get_face_flux();

//                std::cout<<"Face flux in div_u function is:"<<face_flux<<std::endl;
//                sum_face_flux = sum_face_flux + face_flux;
//             }
//   div_u[i] = sum_face_flux;
            
//       }

 for(int i = 0; i<list_of_faces.size(); i++)
       {   
          int owner_index = list_of_faces[i].get_owner_index();
          int neighbour_index = list_of_faces[i].get_neighbour_index();         
          
          double contribution_to_owner = list_of_faces[i].get_face_flux();
          double conttribution_to_neighbour = -1*contribution_to_owner;
          div_u[owner_index] += contribution_to_owner;

          if(neighbour_index > 0)
         {
          div_u[neighbour_index] += conttribution_to_neighbour;
         }

            
       }

       for(int j=0; j<list_of_cells.size(); j++)
       {
         //std::cout<<"div_u_new"<<div_u[j]<<std::endl;

         if(div_u[j] > 1e-6)
         {
          std::cout<<"Alert!!!"<<std::endl;
         }
      }
 } 


void PV_coupling::velocity_output_vector_matrix_coefficients_to_file(double total_cells)
{
     std::ofstream convection_coeffs_a("vel_convection_coeffs_a.txt");
     std::ofstream convection_coeffs_b_x("vel_convection_coeffs_b_x.txt");
     std::ofstream convection_coeffs_b_y("vel_convection_coeffs_b_y.txt");
     std::ofstream diffusion_coeffs_a("vel_diffusion_coeffs_a.txt");
     std::ofstream diffusion_coeffs_b_x("vel_diffusion_coeffs_b_x.txt");
     std::ofstream diffusion_coeffs_b_y("vel_diffusion_coeffs_b_y.txt");
     std::ofstream roc_coeffs_a("vel_roc_coeffs_a.txt");
     std::ofstream roc_coeffs_b_x("vel_roc_coeffs_b_x.txt");
     std::ofstream roc_coeffs_b_y("vel_roc_coeffs_b_y.txt");
     std::ofstream source_coeffs_b_x("vel_source_coeffs_b_x.txt");
     std::ofstream source_coeffs_b_y("vel_source_coeffs_b_y.txt");
     std::ofstream combined_coeffs_a("vel_combined_coeffs_a.txt");
     std::ofstream combined_coeffs_b_x("vel_combined_coeffs_b_x.txt");
     std::ofstream combined_coeffs_b_y("vel_combined_coeffs_b_y.txt");
    
     for(int i = 0; i<total_cells; i++)
     {
        for(int j=0; j<total_cells; j++)
              {
                if(j < total_cells - 1)
                {
                 convection_coeffs_a<<velocity_a_matrix_convection(i,j)<<"\t";
                 diffusion_coeffs_a<<velocity_a_matrix_diffusion(i,j)<<"\t";
                 roc_coeffs_a<<velocity_a_matrix_rate_of_change(i,j)<<"\t";
                 combined_coeffs_a<<velocity_a_matrix_combined(i,j)<<"\t";

                }

                else
                {
                 convection_coeffs_a<<velocity_a_matrix_convection(i,j)<<"\n";
                 diffusion_coeffs_a<<velocity_a_matrix_diffusion(i,j)<<"\n";
                 roc_coeffs_a<<velocity_a_matrix_rate_of_change(i,j)<<"\n";
                 combined_coeffs_a<<velocity_a_matrix_combined(i,j)<<"\n";

                }
              }
     }

          for(int k=0; k < total_cells; k++)
          {
                 convection_coeffs_b_x<<velocity_b_vector_convection_ux(k)<<"\n";
                 convection_coeffs_b_y<<velocity_b_vector_convection_uy(k)<<"\n";

                 diffusion_coeffs_b_x<<velocity_b_vector_diffusion_ux(k)<<"\n";
                 diffusion_coeffs_b_y<<velocity_b_vector_diffusion_uy(k)<<"\n";

                 roc_coeffs_b_x<<velocity_b_vector_rate_of_change_ux(k)<<"\n";
                 roc_coeffs_b_y<<velocity_b_vector_rate_of_change_uy(k)<<"\n";

                 source_coeffs_b_x<<velocity_b_vector_source_ux(k)<<"\n";
                 source_coeffs_b_y<<velocity_b_vector_source_ux(k)<<"\n";

                 combined_coeffs_b_x<<velocity_b_vector_combined_ux(k)<<"\n";
                 combined_coeffs_b_y<<velocity_b_vector_combined_uy(k)<<"\n";
          }    

   convection_coeffs_a.close();
   convection_coeffs_b_x.close();
   convection_coeffs_b_y.close();
   diffusion_coeffs_a.close();
   diffusion_coeffs_b_x.close();
   diffusion_coeffs_b_y.close();
   roc_coeffs_a.close();
   roc_coeffs_b_x.close();
   roc_coeffs_b_y.close();
   combined_coeffs_a.close();
   combined_coeffs_b_x.close();
   combined_coeffs_b_y.close();
   source_coeffs_b_x.close();
   source_coeffs_b_y.close();
}


void PV_coupling::velocity_output_vector_field_to_file(std::vector<double> x_distance, std::vector<double> y_distance)
{
     std::ofstream vector_field_profiles_x("velocity_field_x.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {  
                vector_field_profiles_x<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<velocity.vector_field_values_x(j + i*x_distance.size())<<std::endl;
              }
         }

        std::ofstream vector_field_profiles_y("velocity_field_y.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                vector_field_profiles_y<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<velocity.vector_field_values_y(j + i*x_distance.size())<<std::endl;
              }
         }

        std::ofstream vector_field_profiles_xy("velocity_magnitude.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                 double vx = velocity.vector_field_values_x(j + i*x_distance.size());
                 double vy = velocity.vector_field_values_y(j + i*x_distance.size());
                 double v_mag = sqrt(vx*vx +  vy*vy);

                vector_field_profiles_xy<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<v_mag<<std::endl;
              }
         }


        std::ofstream vx_horizontal("vx_horizontal.txt");

        for(int j=0; j<x_distance.size(); j++)
         {
            vx_horizontal<<x_distance[j]<<"\t"<< velocity.vector_field_values_x(((y_distance.size()/2) - 1.0)*x_distance.size() + j)<<std::endl;
         }


        std::ofstream vx_vertical("vx_vertical.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
            vx_vertical<<velocity.vector_field_values_x((x_distance.size()/2) + i*x_distance.size())<<"\t"<<y_distance[i]<<std::endl;
         }


        std::ofstream vy_horizontal("vy_horizontal.txt");

        for(int j=0; j<x_distance.size(); j++)
         {
            vy_horizontal<<x_distance[j] <<"\t" <<velocity.vector_field_values_y(((y_distance.size()/2) - 1.0)*x_distance.size() + j)<<std::endl;
         }

        std::ofstream vy_vertical("vy_vertical.txt") ;

        for(int i=0; i<y_distance.size(); i++)
         {
            vy_vertical<<y_distance[i]<<"\t"<< velocity.vector_field_values_y((x_distance.size()/2) + i*x_distance.size())<<std::endl;
         }

         std::ofstream vector_plots("vector_plots.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                 double vx = velocity.vector_field_values_x(j + i*x_distance.size());
                 double vy = velocity.vector_field_values_y(j + i*x_distance.size());
                 double v_mag = sqrt(vx*vx +  vy*vy);

                vector_plots<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<pressure.scalar_field_values(j + i*x_distance.size())<<"\t"<<vx/v_mag<<"\t"<<vy/v_mag<<std::endl;
              }
         }


        vector_field_profiles_x.close();
        vector_field_profiles_y.close();
        vector_field_profiles_xy.close();
        vector_plots.close();
        vx_horizontal.close();
        vx_vertical.close();
        vy_horizontal.close();
        vy_vertical.close();
} 

void PV_coupling::velocity_plot_convergence_initial_x(std::ofstream &x_velocity_convergence_initial, int iteration_no)
{   
    int index;
    Eigen::VectorXd diff = (velocity.vector_field_values_x - velocity.vector_old_field_values_x);

    for(int i=0; i<velocity.vector_field_values_x.size(); i++)
    {
        diff(i) = abs(diff(i));
    }

     double max = diff(0);
    
    for(int i=1; i<velocity.vector_field_values_x.size(); i++)
    {   double next_term = diff(i);
        double last_term = diff(i-1);
        
        if(next_term > last_term)
        {
         max = diff(i);
         index = i;
        }
    }

     x_velocity_convergence_initial<<iteration_no<<"\t"<<max<<"\t"<<std::endl;
}

void PV_coupling::velocity_plot_convergence_initial_y(std::ofstream &y_velocity_convergence_initial, int iteration_no)
{
    Eigen::VectorXd difference = (velocity.vector_field_values_y - velocity.vector_old_field_values_y);

    for(int i=0; i<velocity.vector_field_values_y.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<velocity.vector_field_values_y.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }

     y_velocity_convergence_initial<<iteration_no<<"\t"<<max<<std::endl;

}


void PV_coupling::velocity_plot_convergence_final_x(std::ofstream &x_velocity_convergence_final, int iteration_no)
{
       Eigen::VectorXd difference = (velocity.vector_field_values_x - velocity.vector_old_field_values_x);

    for(int i=0; i<velocity.vector_field_values_x.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<velocity.vector_field_values_x.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }

     x_velocity_convergence_final<<iteration_no<<"\t"<<max<<std::endl;
}


void PV_coupling::velocity_plot_convergence_final_y(std::ofstream &y_velocity_convergence_final, int iteration_no)
{
      Eigen::VectorXd difference = (velocity.vector_field_values_y - velocity.vector_old_field_values_y);

    for(int i=0; i<velocity.vector_field_values_y.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<velocity.vector_field_values_y.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }

     y_velocity_convergence_final<<iteration_no<<"\t"<<max<<std::endl;
}


void PV_coupling::pressure_reset_rate_of_change()
{
    pressure_a_matrix_rate_of_change.setZero();
    pressure_b_vector_rate_of_change.setZero();
}

void PV_coupling::pressure_reset_source()
{
    pressure_b_vector_source.setZero();
}

void PV_coupling::pressure_reset_diffusion()
{
   pressure_a_matrix_diffusion.setZero();
   pressure_b_vector_diffusion.setZero();
}

void PV_coupling::pressure_reset_convection()
{
    pressure_a_matrix_convection.setZero();
    pressure_b_vector_convection.setZero();
}

void PV_coupling::pressure_reset_combined()
{
    pressure_a_matrix_combined.setZero();
    pressure_b_vector_combined.setZero();
}

void PV_coupling::pressure_update_a_matrix_rate_of_change(double a, int index)
{
   pressure_a_matrix_rate_of_change(index, index) = pressure_a_matrix_rate_of_change(index, index) + a ;         
}

void PV_coupling::pressure_update_b_vector_rate_of_change(double a, int index)
{
    pressure_b_vector_rate_of_change(index) = pressure_b_vector_rate_of_change(index) + a;
}

void PV_coupling::pressure_update_b_vector_source(double a, int index)
{
    pressure_b_vector_source(index) = pressure_b_vector_source(index) + a;
}

void PV_coupling::pressure_update_a_matrix_diffusion(double a, int row_index, int column_index)
{ 
    pressure_a_matrix_diffusion(row_index, column_index) = pressure_a_matrix_diffusion(row_index, column_index) + a;                
}                                                                         

void PV_coupling::pressure_update_b_vector_diffusion(double a, int index)
{
    pressure_b_vector_diffusion(index) = a;           
}

void PV_coupling::pressure_update_a_matrix_convection(double a, int row_index, int column_index)
{ 
    pressure_a_matrix_convection(row_index, column_index) = pressure_a_matrix_convection(row_index, column_index) + a;                 
}                                                                         

void PV_coupling::pressure_update_b_vector_convection(double a, int index)
{
    pressure_b_vector_convection(index) = pressure_b_vector_convection(index) + a;           
}

void PV_coupling::pressure_combine_a_and_b_matrices()
{
    for(int i = 0; i<pressure_a_matrix_combined.rows(); i++)
    {
        for(int j=0; j<pressure_a_matrix_combined.cols(); j++)
        {
            pressure_a_matrix_combined(i,j) = pressure_a_matrix_rate_of_change(i,j) + pressure_a_matrix_diffusion(i,j) + pressure_a_matrix_convection(i,j);
        }
    }

    for(int k =0 ; k <pressure_b_vector_combined.size(); k++)
    {
        pressure_b_vector_combined(k) = pressure_b_vector_rate_of_change(k) + pressure_b_vector_source(k) + pressure_b_vector_diffusion(k) + pressure_b_vector_convection(k) ;
    }
}

void PV_coupling::pressure_compute_source_matrix(std::vector<Cell> list_of_cells, std::vector<Face> list_of_faces)                                                               
 {  
    pressure_reset_source(); 

    for(int i = 0; i<list_of_faces.size(); i++)
    {  
        int face_index = list_of_faces[i].get_face_index();
        int owner_index = list_of_faces[i].get_owner_index();
        int neighbour_index = list_of_faces[i].get_neighbour_index();

        double face_flux_value = list_of_faces[i].get_face_flux();

        double face_flux_value_owner_cell = face_flux_value;
        double face_flux_value_neighbour_cell = -1 * face_flux_value;

        if(list_of_faces[i].get_is_boundary_face() == true)
        { 
          pressure_update_b_vector_source(face_flux_value_owner_cell, owner_index);
        }

        else
        {
           pressure_update_b_vector_source(face_flux_value_owner_cell, owner_index);
           pressure_update_b_vector_source(face_flux_value_neighbour_cell, neighbour_index);
        } 
    }
 } 

void PV_coupling::pressure_compute_diffusion_matrix(std::vector<Face> list_of_faces,  std::vector <double> ap_coeff, std::vector<Boundary> list_of_boundaries, int iteration_no, Scalar_boundary_field pressure_boundary)                                                              
 { 
      pressure_reset_diffusion();  

      for(int i = 0; i<list_of_faces.size(); i++)
      {
          double a_owner_diagonal, a_neighbour_diagonal, a_owner_neighbour, a_neighbour_owner, b_owner_source, b_neighbour_source;
            
            int face_owner_index = list_of_faces[i].get_owner_index();
            int face_neighbour_index = list_of_faces[i].get_neighbour_index();

            double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
            double face_delta = list_of_faces[i].get_face_delta(); 

            Vector face_normal = list_of_faces[i].get_face_normal();

            double face_normal_magnitude = face_normal.magnitude_of_vector();

                if(list_of_faces[i].get_is_boundary_face() == true)
                {
                     int boundary_index = list_of_faces[i].get_boundary_index();
                     Boundary b = list_of_boundaries[boundary_index];

                     //if(b.return_boundary_name() == "empty")
                     if(pressure_boundary.get_scalar_boundary_name(boundary_index) == "empty")
                     {
                        continue;
                     }

                     else
                     {
                            double ap_coeff_owner = ap_coeff[face_owner_index];
                            double ap_boundary_face = interpolation_factor * ap_coeff_owner;

                        // if(b.get_is_scalar_dirichlet() == true)
                        if(pressure_boundary.get_scalar_boundary_type(boundary_index) == "dirichlet")
                        {   
                            std::vector<double> dirichlet_value = pressure_boundary.get_scalar_boundary_field_values();

                            a_owner_diagonal = -1*ap_boundary_face * face_normal_magnitude * face_delta;
                            pressure_update_a_matrix_diffusion(a_owner_diagonal, face_owner_index, face_owner_index); 


                            b_owner_source = -1* ap_boundary_face * face_normal_magnitude * face_delta * dirichlet_value[boundary_index];
                            pressure_update_b_vector_diffusion(b_owner_source, face_owner_index);
                        }

                        else
                        {
                            std::vector<double> neumann_gradient =  pressure_boundary.get_scalar_boundary_field_values();

                            b_owner_source = -1 * ap_coeff_owner * face_normal_magnitude * neumann_gradient[boundary_index] ;

                            pressure_update_b_vector_diffusion(b_owner_source, face_owner_index);
                        }
                     }
                }

                else
                {
                   double ap_owner_cell = ap_coeff[face_owner_index];
                   double ap_neighbour_cell = ap_coeff[face_neighbour_index];

                   double ap_face = (interpolation_factor*ap_owner_cell) + ((1 - interpolation_factor)*(ap_neighbour_cell)) ;

                   a_owner_diagonal = -1 * ap_face * face_normal_magnitude * face_delta;
                   a_neighbour_diagonal = -1 * ap_face * face_normal_magnitude * face_delta;

                   a_owner_neighbour =  1 * ap_face * face_normal_magnitude * face_delta;
                   a_neighbour_owner =  1 * ap_face * face_normal_magnitude * face_delta;
                   
                   pressure_update_a_matrix_diffusion(a_owner_diagonal, face_owner_index, face_owner_index);
                   pressure_update_a_matrix_diffusion(a_neighbour_diagonal, face_neighbour_index, face_neighbour_index);
                   pressure_update_a_matrix_diffusion(a_owner_neighbour, face_owner_index, face_neighbour_index);
                   pressure_update_a_matrix_diffusion(a_neighbour_owner, face_neighbour_index, face_owner_index);
            
                }
      }
 } 
    

void PV_coupling::pressure_calculate_initial_residuals_p(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_initial_residual_p)
{
     if(n == 199)
  { 
   std::ofstream file_initial_residuals_p("pressure_initial_residuals.txt");

   pressure_initial_residuals_p = pressure_b_vector_combined - pressure_a_matrix_combined * pressure.scalar_field_values;
   
   for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
         file_initial_residuals_p<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<pressure.scalar_field_values(j + (i*x_distance.size()))<<std::endl;
            
      }
   } 

   file_initial_residuals_p.close();
  } 
 
  pressure_initial_residuals_p = pressure_b_vector_combined - pressure_a_matrix_combined * pressure.scalar_field_values;

  double sum_r_p = 0.0;

  for(int i=0; i<pressure_initial_residuals_p.size(); i++)
  {
     sum_r_p = sum_r_p + pressure_initial_residuals_p(i);
  } 

  sum_initial_residual_p<<n<<"\t"<<sum_r_p<<std::endl;

}
        
void PV_coupling::pressure_combine_and_solve_matrices(std::vector<Cell> list_of_cells)
{
    pressure_combine_a_and_b_matrices();
    int size = list_of_cells.size();

    pressure.scalar_field_values = pressure_a_matrix_combined.fullPivHouseholderQr().solve(pressure_b_vector_combined);

    pressure.set_old_scalar_field_values();

}


void PV_coupling::pressure_under_relax()
{
    pressure.scalar_field_values = pressure.scalar_old_field_values + (0.7) * (pressure.scalar_field_values - pressure.scalar_old_field_values); 
} 

void PV_coupling::pressure_calculate_final_residuals_p(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_final_residual_p)
{
  if(n == 199)
  { 
   std::ofstream file_final_residuals_p("pressure_final_residuals.txt");

   pressure_final_residuals_p = pressure_b_vector_combined - pressure_a_matrix_combined * pressure.scalar_field_values;
   
   for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
            file_final_residuals_p<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<pressure.scalar_field_values(j + (i*x_distance.size()))<<std::endl;     
      }
   } 

   file_final_residuals_p.close();
  }


  pressure_final_residuals_p = pressure_b_vector_combined - pressure_a_matrix_combined * pressure.scalar_field_values;

  double sum_r_p = 0.0;

  for(int i=0; i<pressure_final_residuals_p.size(); i++)
  {
     sum_r_p = sum_r_p + pressure_final_residuals_p(i);
  } 

  sum_final_residual_p<<n<<"\t"<<sum_r_p<<std::endl;
}


// void PV_coupling::pressure_compute_flux_correction(std::vector<Cell> &list_of_cells, std::vector<Face> &list_of_faces, std::vector<double> temp_ap_coeffs, std::vector<Boundary> list_of_boundaries)
// {       
//    double ap_face, grad_p_face, new_face_flux_value;

//    for(int i=0; i<list_of_faces.size(); i++)
//    {
//       int owner_index = list_of_faces[i].get_owner_index();
//       int neighbour_index = list_of_faces[i].get_neighbour_index();
//       int interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
//       int face_delta = list_of_faces[i].get_face_delta();

//       Vector face_normal = list_of_faces[i].get_face_normal();

//       double face_normal_magnitude = face_normal.magnitude_of_vector();

//       if(list_of_faces[i].get_is_boundary_face() == true)
//         {   
//             int boundary_index = list_of_faces[i].get_boundary_index();
//             Boundary b = list_of_boundaries[boundary_index];

//             if(b.return_boundary_name() == "empty")
//             {
//             continue;
//             }

//             else
//             {
//                 double IF = list_of_faces[i].get_face_interpolation_factor();
//                 double ap_owner = temp_ap_coeffs[owner_index];

//                 grad_p_face = 0;
//                 new_face_flux_value = ap_face * grad_p_face;
//             }
//         }

//         else
//         {
//             double ap_owner_cell = temp_ap_coeffs[owner_index];
//             double ap_neighbour_cell = temp_ap_coeffs[neighbour_index];

//             ap_face = (interpolation_factor*ap_owner_cell) + ((1 - interpolation_factor)*(ap_neighbour_cell)) ;

//             grad_p_face = ((pressure.scalar_field_values(neighbour_index) - pressure.scalar_field_values(owner_index))*face_delta); 

//             new_face_flux_value =  ap_face* face_normal_magnitude *grad_p_face;
            
//         }
//               list_of_faces[i].set_face_flux(new_face_flux_value); 
//    }
// }

void PV_coupling::pressure_compute_flux_correction(std::vector<Cell> &list_of_cells, std::vector<Face> &list_of_faces, std::vector<double> temp_ap_coeffs, std::vector<Boundary> list_of_boundaries)
{        
   double ap_face, grad_p_face, new_face_flux_value, f_corr, test_val = 0.0;

   for(int i=0; i<list_of_faces.size(); i++)
   {
      int owner_index = list_of_faces[i].get_owner_index();
      int neighbour_index = list_of_faces[i].get_neighbour_index();
      double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
      double face_delta = list_of_faces[i].get_face_delta();

      Vector face_normal = list_of_faces[i].get_face_normal();

      double face_normal_magnitude = face_normal.magnitude_of_vector();

      if(list_of_faces[i].get_is_boundary_face() == true)
        {   
            int boundary_index = list_of_faces[i].get_boundary_index();
            Boundary b = list_of_boundaries[boundary_index];

            if(b.return_boundary_name() == "empty")
            {
            continue;
            }

            else
            {
                double IF = list_of_faces[i].get_face_interpolation_factor();
                double ap_owner = temp_ap_coeffs[owner_index];

                grad_p_face = 0;
                new_face_flux_value = ap_face * grad_p_face;

                f_corr = list_of_faces[i].get_face_flux() - new_face_flux_value;
            }
        }

        else
        {
            double ap_owner_cell = temp_ap_coeffs[owner_index];
            double ap_neighbour_cell = temp_ap_coeffs[neighbour_index];

            ap_face = (interpolation_factor*ap_owner_cell) + ((1 - interpolation_factor)*(ap_neighbour_cell)) ;

            grad_p_face = ((pressure.scalar_field_values(neighbour_index) - pressure.scalar_field_values(owner_index))*face_delta); 

            new_face_flux_value =  ap_face* face_normal_magnitude *grad_p_face;
             f_corr = list_of_faces[i].get_face_flux() - new_face_flux_value;
        }
            list_of_faces[i].set_face_flux(f_corr); 

   }
}




          
void PV_coupling::pressure_output_scalar_matrix_coefficients_to_file(double total_cells)
{
     std::ofstream pressure_convection_coeffs_a("pressure_convection_coeffs_a.txt");
     std::ofstream pressure_convection_coeffs_b("pressure_convection_coeffs_b.txt");
     std::ofstream pressure_diffusion_coeffs_a("pressure_diffusion_coeffs_a.txt");
     std::ofstream pressure_diffusion_coeffs_b("pressure_diffusion_coeffs_b.txt");
     std::ofstream pressure_roc_coeffs_a("pressure_roc_coeffs_a.txt");
     std::ofstream pressure_roc_coeffs_b("pressure_roc_coeffs_b.txt");
     std::ofstream pressure_source_coeffs_b("pressure_source_coeffs_b.txt");
     std::ofstream pressure_combined_coeffs_a("pressure_combined_coeffs_a.txt");
     std::ofstream pressure_combined_coeffs_b("pressure_combined_coeffs_b.txt");
    
     for(int i = 0; i<total_cells; i++)
     {
        for(int j=0; j<total_cells; j++)
              {
                if(j < total_cells - 1)
                {
                 pressure_convection_coeffs_a<<pressure_a_matrix_convection(i,j)<<"\t";
                 pressure_diffusion_coeffs_a<<pressure_a_matrix_diffusion(i,j)<<"\t";
                 pressure_roc_coeffs_a<<pressure_a_matrix_rate_of_change(i,j)<<"\t";
                 pressure_combined_coeffs_a<<pressure_a_matrix_combined(i,j)<<"\t";

                }

                else
                {
                 pressure_convection_coeffs_a<<pressure_a_matrix_convection(i,j)<<"\n";
                 pressure_diffusion_coeffs_a<<pressure_a_matrix_diffusion(i,j)<<"\n";
                 pressure_roc_coeffs_a<<pressure_a_matrix_rate_of_change(i,j)<<"\n";
                 pressure_combined_coeffs_a<<pressure_a_matrix_combined(i,j)<<"\n";

                }
              }
     }

     for(int k=0; k<total_cells; k++)
     {
                pressure_convection_coeffs_b<<pressure_b_vector_convection(k)<<"\n";
                 pressure_diffusion_coeffs_b<<pressure_b_vector_diffusion(k)<<"\n";
                 pressure_roc_coeffs_b<<pressure_b_vector_rate_of_change(k)<<"\n";
                 pressure_source_coeffs_b<<pressure_b_vector_source(k)<<"\n";
                 pressure_combined_coeffs_b<<pressure_b_vector_combined(k)<<"\n";

     }

   pressure_convection_coeffs_a.close();
   pressure_convection_coeffs_b.close();
   pressure_diffusion_coeffs_a.close();
   pressure_diffusion_coeffs_b.close();
   pressure_roc_coeffs_a.close();
   pressure_roc_coeffs_b.close();
   pressure_combined_coeffs_a.close();
   pressure_combined_coeffs_b.close();
   pressure_source_coeffs_b.close();
}


void PV_coupling::pressure_output_scalar_field_to_file(std::vector<double> x_distance, std::vector<double> y_distance, std::vector<Cell> list_of_cells)
{
     std::ofstream scalar_field_profiles("pressure.txt");
     std::ofstream scalar_field_profiles_volume("Pressure_vol.txt");
     std::ofstream pressure_x_line("pressure_x_line.txt");
     std::ofstream pressure_y_line("pressure_y_line.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                scalar_field_profiles_volume<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<pressure.scalar_field_values(j + i*x_distance.size())/cell_v<<std::endl;
                scalar_field_profiles<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<pressure.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }

        for(int j=0; j<x_distance.size(); j++)
         {
            pressure_x_line<<x_distance[j]<<"\t"<< pressure.scalar_field_values(((y_distance.size()/2) - 1.0)*x_distance.size() + j)<<std::endl;
         }

         for(int i=0; i<y_distance.size(); i++)
         {
            pressure_y_line<<y_distance[i]<<"\t"<<pressure.scalar_field_values((x_distance.size()/2) + i*x_distance.size())<<std::endl;
         }

        scalar_field_profiles.close();
        pressure_x_line.close();
        pressure_y_line.close();
}     
                
 void PV_coupling::pressure_plot_convergence_initial(std::ofstream &plot_convergence_initial_h, int iteration_no)
 {
    Eigen::VectorXd difference = (pressure.scalar_field_values - pressure.scalar_old_field_values);

    for(int i=0; i<pressure.scalar_field_values.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<pressure.scalar_field_values.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }
     plot_convergence_initial_h<<iteration_no<<"\t"<<max<<std::endl;
 }               



void PV_coupling::pressure_compute_source_matrix_vol_based(std::vector<Cell> list_of_cells)                                                               
 {  
    pressure_reset_source(); 
    
    for(int i = 0; i<list_of_cells.size(); i++)
    {  
        int cell_index = list_of_cells[i].get_cell_index();

        double source_value = list_of_cells[i].get_sum_of_fluxes_through_cell();

        pressure_update_b_vector_source(source_value, cell_index);
    }
 } 


Scalar_field PV_coupling::retrieve_pressure_field()                                                               
 {  
     return pressure;
 } 


 void PV_coupling::set_alpha_scalar_initial_fields(std::vector<Cell> list_of_cells, double rho_1, double rho_2, double nu_1, double nu_2, double mu_1, double mu_2, double x_max, double y_max)
{
    std::ofstream plot_initial_alpha("initial_alpha.txt");
    std::ofstream plot_initial_nu("initial_nu.txt");
    std::ofstream plot_initial_rho("initial_rho.txt");
    std::ofstream plot_initial_mu("initial_mu.txt");

    for(int i=0; i<list_of_cells.size(); i++)
    {
        Vector distance = list_of_cells[i].get_cell_centre();

        // if(distance[0]<= 0.025 && distance[1] <= 0.05)
        if((distance[0]<= 0.075 && distance[0] >= 0.025)  && (distance[1] <= 0.075 && distance[1] >= 0.025))
        {
            alpha.scalar_field_values(i) = 1.0;
        }

        else
        {
            alpha.scalar_field_values(i) = 0.0;
        }

        

        rho.scalar_field_values(i) = alpha.scalar_field_values(i) * rho_1 + ((1 - (alpha.scalar_field_values(i))) * rho_2);
        nu.scalar_field_values(i) = alpha.scalar_field_values(i) * nu_1 + ((1 - alpha.scalar_field_values(i)) * nu_2);
        mu.scalar_field_values(i) = alpha.scalar_field_values(i) * mu_1 + ((1 - alpha.scalar_field_values(i)) * mu_2);

        plot_initial_alpha<<distance[0]<<"\t"<<distance[1]<<"\t"<<alpha.scalar_field_values(i)<<std::endl;
        plot_initial_nu<<distance[0]<<"\t"<<distance[1]<<"\t"<<nu.scalar_field_values(i)<<std::endl;
        plot_initial_rho<<distance[0]<<"\t"<<distance[1]<<"\t"<<rho.scalar_field_values(i)<<std::endl;
        plot_initial_mu<<distance[0]<<"\t"<<distance[1]<<"\t"<<mu.scalar_field_values(i)<<std::endl;
    }

    plot_initial_alpha.close();
    plot_initial_nu.close();
    plot_initial_rho.close();
    plot_initial_mu.close();
}

void PV_coupling::alpha_reset_rate_of_change()
{
    alpha_a_matrix_rate_of_change.setZero();
    alpha_b_vector_rate_of_change.setZero();
}

void PV_coupling::alpha_reset_source()
{
    alpha_b_vector_source.setZero();
}

void PV_coupling::alpha_reset_convection()
{
    alpha_a_matrix_convection.setZero();
    alpha_b_vector_convection.setZero();
}

void PV_coupling::alpha_reset_combined()
{
    alpha_a_matrix_combined.setZero();
    alpha_b_vector_combined.setZero();
}

void PV_coupling::alpha_update_a_matrix_rate_of_change(double a, int index)
{
   alpha_a_matrix_rate_of_change(index, index) = alpha_a_matrix_rate_of_change(index, index) + a ;         
}

void PV_coupling::alpha_update_b_vector_rate_of_change(double a, int index)
{
    alpha_b_vector_rate_of_change(index) = alpha_b_vector_rate_of_change(index) + a;
}

void PV_coupling::alpha_update_b_vector_source(double a, int index)
{
    alpha_b_vector_source(index) = alpha_b_vector_source(index) + a;
}

void PV_coupling::alpha_update_a_matrix_convection(double a, int row_index, int column_index)
{ 
    alpha_a_matrix_convection(row_index, column_index) = alpha_a_matrix_convection(row_index, column_index) + a;                 
}                                                                         

void PV_coupling::alpha_update_b_vector_convection(double a, int index)
{
    alpha_b_vector_convection(index) = alpha_b_vector_convection(index) + a;           
}

void PV_coupling::alpha_combine_a_and_b_matrices()
{
    for(int i = 0; i<alpha_a_matrix_combined.rows(); i++)
    {
        for(int j=0; j<alpha_a_matrix_combined.cols(); j++)
        {
            alpha_a_matrix_combined(i,j) = alpha_a_matrix_rate_of_change(i,j) + alpha_a_matrix_convection(i,j);
        }
    }

    for(int k =0 ; k <alpha_b_vector_combined.size(); k++)
    {
        alpha_b_vector_combined(k) = alpha_b_vector_rate_of_change(k) + alpha_b_vector_source(k) + alpha_b_vector_convection(k) ;
    }
}

 void PV_coupling::alpha_rate_of_change_discretization(std::vector<Cell> list_of_cells, double delta_time)                                                               
 {
      //std::cout<<"Running";
      alpha_reset_rate_of_change();  

      for(int i = 0; i<list_of_cells.size(); i++)
      {  
            int cell_index = list_of_cells[i].get_cell_index();
            double cell_vol = list_of_cells[cell_index].get_cell_volume();
            double alpha_old = alpha.scalar_field_values(cell_index);

            double rate_of_change_diagonal_contribution = cell_vol/delta_time;
            double rate_of_change_source_contribution = (cell_vol * alpha_old)/ delta_time;

            alpha_update_a_matrix_rate_of_change(rate_of_change_diagonal_contribution, cell_index);
            alpha_update_b_vector_rate_of_change(rate_of_change_source_contribution, cell_index);
      }
 } 

 void PV_coupling::alpha_convection_discretization(std::vector<Face> list_of_faces, std::vector<Cell> list_of_cells, std::vector<Boundary> list_of_boundaries, Scalar_boundary_field alpha_boundary, Vector_boundary_field velocity_boundary)                                                               
 { 
      alpha_reset_convection();  

      Vector_field grad_alpha = alpha.compute_gauss_gradient(list_of_cells, list_of_faces);

      for(int i = 0; i<list_of_faces.size(); i++)
        {   
             double face_alpha, velocity_flux, beta = 0.5, face_area, phi_c_tilda;
             Vector face_velocity, face_normal; 

             int face_index = list_of_faces[i].get_face_index();

             double a_owner_diagonal, a_neighbour_diagonal, a_owner_neighbour, a_neighbour_owner, b_owner_source, b_neighbour_source;
             double neighbour_velocity_x, neighbour_velocity_y, neighbour_velocity_z;

             int owner_index = list_of_faces[i].get_owner_index();
             int neighbour_index = list_of_faces[i].get_neighbour_index();

             double interpolation_factor = list_of_faces[i].get_face_interpolation_factor();
             double face_delta = list_of_faces[i].get_face_interpolation_factor();
          
             face_normal = list_of_faces[i].get_face_normal();
             face_area = face_normal.magnitude_of_vector();   

            if(list_of_faces[i].get_is_boundary_face() == true)
            {
                int boundary_index = list_of_faces[i].get_boundary_index();
                Boundary b = list_of_boundaries[list_of_faces[i].get_boundary_index()];

                std::vector<std::vector<double>> vel = velocity_boundary.get_vector_boundary_field_values();

                face_velocity = Vector(vel[boundary_index][0], vel[boundary_index][1], vel[boundary_index][2]);
            }

            else
            {      
                double owner_vx = velocity.vector_field_values_x(owner_index);
                double owner_vy = velocity.vector_field_values_y(owner_index);
                double owner_vz = velocity.vector_field_values_z(owner_index);
                double neighbour_vx = velocity.vector_field_values_x(neighbour_index);
                double neighbour_vy = velocity.vector_field_values_y(neighbour_index);
                double neighbour_vz = velocity.vector_field_values_z(neighbour_index);
                Vector owner_velocity(owner_vx, owner_vy, owner_vz);
                Vector neighbour_velocity(neighbour_vx, neighbour_vy, neighbour_vz);
                face_velocity = owner_velocity * interpolation_factor + neighbour_velocity * (1 - interpolation_factor);

                int base_cell_index = owner_index;
                int ud_cell_index = neighbour_index;

                double alpha_base_cell = alpha.scalar_field_values(base_cell_index);
                double alpha_ud_cell = alpha.scalar_field_values(ud_cell_index);
                double alpha_diff = alpha_ud_cell - alpha_base_cell;

                //std::cout<<"alpha_diff: "<<alpha_diff<<std::endl;

                Vector alpha_owner_cell_gradient = Vector(grad_alpha.vector_field_values_x(base_cell_index), grad_alpha.vector_field_values_y(base_cell_index), grad_alpha.vector_field_values_z(base_cell_index));
                //std::cout<<"alpha_ownwe_grad:"<<alpha_owner_cell_gradient[0]<<"y:  "<<alpha_owner_cell_gradient[1]<<"z:" <<alpha_owner_cell_gradient[2]<<std::endl;

                Vector base_cell_centre = list_of_cells[base_cell_index].get_cell_centre();
                Vector ud_cell_centre = list_of_cells[ud_cell_index].get_cell_centre();

                Vector cell_to_cell_distance = base_cell_centre - ud_cell_centre ;

                double grad_phi_dot_p = alpha_owner_cell_gradient.dot_product(cell_to_cell_distance);
                double phi_tilda_denominator = 2.0*grad_phi_dot_p;  
                //std::cout<<"phi_t_deno: "<<phi_tilda_denominator<<std::endl;
         
                if(grad_phi_dot_p < 1e-30)
                {
                  phi_c_tilda = 1;
                }

                else
                {
                   phi_c_tilda = 1.0 - (alpha_diff/phi_tilda_denominator);
                }

                
                //std::cout<<phi_c_tilda<<std::endl;
            }


               velocity_flux = face_normal.dot_product(face_velocity);

            if((phi_c_tilda > 0) && (phi_c_tilda < 0.5))
            {
              double gamma = (phi_c_tilda / beta);
              double interpolation_factor = (1.0 - (gamma/0.5));

              // if(velocity_flux >= 0)
              //  {
              //   interpolation_factor = 1;
              //  }

              //  else
              //  {
              //   interpolation_factor = 0;
              //  }
            } 

            else if((phi_c_tilda >= 0.5) && (phi_c_tilda < 1.0))
            {
               if(velocity_flux >= 0)
               {
                interpolation_factor = 0;
               }

               else
               {
                interpolation_factor = 1;
               }

              // if(velocity_flux >= 0)
              //  {
              //   interpolation_factor = 1;
              //  }

              //  else
              //  {
              //   interpolation_factor = 0;
              //  }
            }

            else
            {
               if(velocity_flux >= 0)
               {
                interpolation_factor = 1;
               }

               else
               {
                interpolation_factor = 0;
               }
            }
          

            if(list_of_faces[i].get_is_boundary_face() == true)
            {
                  
                  int boundary_index = list_of_faces[i].get_boundary_index();
                  Boundary b = list_of_boundaries[boundary_index];

                if(alpha_boundary.get_scalar_boundary_name(boundary_index) == "empty")
                  {
                    continue;
                  }

                  else
                  {

                    if(alpha_boundary.get_scalar_boundary_type(boundary_index) == "dirichlet")
                    {   
                        std::vector<double> dirichlet_value = alpha_boundary.get_scalar_boundary_field_values();

                        b_owner_source = -1*velocity_flux * dirichlet_value[boundary_index];
                        alpha_update_b_vector_convection(b_owner_source, owner_index);
                    }

                    if(alpha_boundary.get_scalar_boundary_type(boundary_index) == "neumann")
                    {
                        std::vector<double> neumann_gradient = alpha_boundary.get_scalar_boundary_field_values();

                        a_owner_diagonal = velocity_flux;
                        alpha_update_a_matrix_convection(a_owner_diagonal, owner_index, owner_index);

                        b_owner_source = -1*(velocity_flux/face_delta) * neumann_gradient[boundary_index];
                        alpha_update_b_vector_convection(b_owner_source, owner_index);
                    }
                  }
            }

            else
            {       

              a_owner_diagonal = interpolation_factor*velocity_flux;
              a_neighbour_diagonal = -1*(1 - interpolation_factor)*velocity_flux;

                alpha_update_a_matrix_convection(a_owner_diagonal, owner_index, owner_index);
                alpha_update_a_matrix_convection(a_neighbour_diagonal, neighbour_index, neighbour_index);

                a_owner_neighbour = (1 - interpolation_factor)*velocity_flux;
                a_neighbour_owner = -1*velocity_flux*interpolation_factor;

              alpha_update_a_matrix_convection(a_owner_neighbour, owner_index, neighbour_index);
              alpha_update_a_matrix_convection(a_neighbour_owner, neighbour_index, owner_index);
                
            } 
        }
    }  


void PV_coupling::alpha_combine_and_solve_matrices()
{
    alpha_combine_a_and_b_matrices();

    alpha.scalar_field_values = alpha_a_matrix_combined.fullPivHouseholderQr().solve(alpha_b_vector_combined);
    //alpha.scalar_field_values = alpha_a_matrix_combined.lu().solve(alpha_b_vector_combined);
}

void PV_coupling::alpha_update_rho_nu(std::vector<Cell> list_of_cells, double rho_1, double rho_2, double nu_1, double nu_2, double mu_1, double mu_2)
{
    rho.set_old_scalar_field_values();
    nu.set_old_scalar_field_values();
    mu.set_old_scalar_field_values();

    for(int i=0; i<list_of_cells.size(); i++)
    {
        int cell_index = list_of_cells[i].get_cell_index();

        double alpha_cell = alpha.scalar_field_values(cell_index);

        rho.scalar_field_values(cell_index) = (alpha_cell * rho_1) + ((1.0 - alpha_cell ) * rho_2);
        nu.scalar_field_values(cell_index) = (alpha_cell * nu_1) + ((1.0 - alpha_cell ) * nu_2);
        mu.scalar_field_values(cell_index) = (alpha_cell * mu_1) + ((1.0 - alpha_cell ) * mu_2);
    }
}

void PV_coupling::alpha_output_scalar_fields_to_file(std::vector<double> x_distance, std::vector<double> y_distance, std::vector<Cell> list_of_cells, int iteration_no)
{
     
    //  std::ofstream alpha_field_profiles_10("alpha_10.txt");
    //  std::ofstream alpha_field_profiles_20("alpha_20.txt");
    //  std::ofstream alpha_field_profiles_30("alpha_30.txt");
    //  std::ofstream alpha_field_profiles_40("alpha_40.txt");
    //  std::ofstream alpha_field_profiles_50("alpha_50.txt");
    //  std::ofstream alpha_field_profiles_60("alpha_60.txt");
    //  std::ofstream alpha_field_profiles_70("alpha_70.txt");
    //  std::ofstream alpha_field_profiles_80("alpha_80.txt");
    //  std::ofstream alpha_field_profiles_90("alpha_90.txt");
    //  std::ofstream alpha_field_profiles_100("alpha_100.txt");
    //  std::ofstream alpha_field_profiles_110("alpha_110.txt");
    //  std::ofstream alpha_field_profiles_120("alpha_120.txt");
    //  std::ofstream alpha_field_profiles_130("alpha_130.txt");
    //  std::ofstream alpha_field_profiles_140("alpha_140.txt");
    //  std::ofstream alpha_field_profiles_150("alpha_150.txt");
    //  std::ofstream alpha_field_profiles_160("alpha_160.txt");
    //  std::ofstream alpha_field_profiles_170("alpha_170.txt");
    //  std::ofstream alpha_field_profiles_180("alpha_180.txt");
    //  std::ofstream alpha_field_profiles_190("alpha_190.txt");
    //  std::ofstream alpha_field_profiles_200("alpha_200.txt");
    //  std::ofstream alpha_field_profiles_210("alpha_210.txt");
    //  std::ofstream alpha_field_profiles_220("alpha_220.txt");
    //  std::ofstream alpha_field_profiles_230("alpha_230.txt");
    //  std::ofstream alpha_field_profiles_240("alpha_240.txt");
    //  std::ofstream alpha_field_profiles_250("alpha_250.txt");
    //  std::ofstream alpha_field_profiles_260("alpha_260.txt");
    //  std::ofstream alpha_field_profiles_270("alpha_270.txt");
    //  std::ofstream alpha_field_profiles_280("alpha_280.txt");
    //  std::ofstream alpha_field_profiles_290("alpha_290.txt");
    //  std::ofstream alpha_field_profiles_300("alpha_300.txt");


    if(iteration_no == 1)
      {  std::ofstream alpha_field_profiles_1("alpha_1.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1.close();
      }  

     if(iteration_no == 10)
      {std::ofstream alpha_field_profiles_10("alpha_10.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_10<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_10.close();
      } 

          if(iteration_no == 20)
      {std::ofstream alpha_field_profiles_20("alpha_20.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_20<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_20.close();
      } 

          if(iteration_no == 30)
      {std::ofstream alpha_field_profiles_30("alpha_30.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_30<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_30.close();
      } 

          if(iteration_no == 40)
      {std::ofstream alpha_field_profiles_40("alpha_40.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_40<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_40.close();
      } 

          if(iteration_no == 50)
      {std::ofstream alpha_field_profiles_50("alpha_50.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_50<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_50.close();
      } 

          if(iteration_no == 60)
      {std::ofstream alpha_field_profiles_60("alpha_60.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_60<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_60.close();
      } 

          if(iteration_no == 70)
      {std::ofstream alpha_field_profiles_70("alpha_70.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_70<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_70.close();
      } 

          if(iteration_no == 80)
      {std::ofstream alpha_field_profiles_80("alpha_80.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_80<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_80.close();
      } 

          if(iteration_no == 90)
      {std::ofstream alpha_field_profiles_90("alpha_90.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_90<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_90.close();
      } 


                if(iteration_no == 100)
      {std::ofstream alpha_field_profiles_100("alpha_100.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_100<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_100.close();
      } 

                if(iteration_no == 110)
      {std::ofstream alpha_field_profiles_110("alpha_110.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_110<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_110.close();
      } 

                if(iteration_no == 120)
      {std::ofstream alpha_field_profiles_120("alpha_120.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_120<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_120.close();
      } 

                if(iteration_no == 130)
      {std::ofstream alpha_field_profiles_130("alpha_130.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_130<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_130.close();
      } 

                if(iteration_no == 140)
      {std::ofstream alpha_field_profiles_140("alpha_140.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_140<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_140.close();
      } 

                if(iteration_no == 150)
      {std::ofstream alpha_field_profiles_150("alpha_150.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_150<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_150.close();
      } 

                if(iteration_no == 160)
      {std::ofstream alpha_field_profiles_160("alpha_160.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_160<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_160.close();
      } 

                if(iteration_no == 170)
      {std::ofstream alpha_field_profiles_170("alpha_170.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_170<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_170.close();
      } 

                if(iteration_no == 180)
      {std::ofstream alpha_field_profiles_180("alpha_180.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_180<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_180.close();
      } 

                if(iteration_no == 190)
      {std::ofstream alpha_field_profiles_190("alpha_190.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_190<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_190.close();
      } 

                if(iteration_no == 200)
      {std::ofstream alpha_field_profiles_200("alpha_200.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_200<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_200.close();
      } 

                      if(iteration_no == 210)
      {std::ofstream alpha_field_profiles_210("alpha_210.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_210<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_210.close();
      } 

                      if(iteration_no == 220)
      {std::ofstream alpha_field_profiles_220("alpha_220.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_220<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_220.close();
      } 

                      if(iteration_no == 230)
      {std::ofstream alpha_field_profiles_230("alpha_230.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_230<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_230.close();
      } 

                      if(iteration_no == 240)
      {std::ofstream alpha_field_profiles_240("alpha_240.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_240<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_240.close();
      } 

                      if(iteration_no == 250)
      {std::ofstream alpha_field_profiles_250("alpha_250.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_250<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_250.close();
      } 

                            if(iteration_no == 260)
      {std::ofstream alpha_field_profiles_260("alpha_260.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_260<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_260.close();
      } 

                            if(iteration_no == 270)
      {std::ofstream alpha_field_profiles_270("alpha_270.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_270<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_270.close();
      } 

                            if(iteration_no == 280)
      {std::ofstream alpha_field_profiles_280("alpha_280.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_280<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_280.close();
      } 

                            if(iteration_no == 290)
      {std::ofstream alpha_field_profiles_290("alpha_290.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_290<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_290.close();
      } 

                            if(iteration_no == 300)
      {std::ofstream alpha_field_profiles_300("alpha_300.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_300<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_300.close();
      } 

                                  if(iteration_no == 400)
      {std::ofstream alpha_field_profiles_400("alpha_400.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_400<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_400.close();
      } 
                                  if(iteration_no == 600)
      {std::ofstream alpha_field_profiles_600("alpha_600.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_600<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_600.close();
      } 
                                  if(iteration_no == 800)
      {std::ofstream alpha_field_profiles_800("alpha_800.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_800<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_800.close();
      } 
                                  if(iteration_no == 1000)
      {std::ofstream alpha_field_profiles_1000("alpha_1000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1000.close();
      } 
                                  if(iteration_no == 1100)
      {std::ofstream alpha_field_profiles_1100("alpha_1100.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1100<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1100.close();
      } 
                                  if(iteration_no == 1200)
      {std::ofstream alpha_field_profiles_1200("alpha_1200.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1200<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1200.close();
      } 
                                  if(iteration_no == 1300)
      {std::ofstream alpha_field_profiles_1300("alpha_1300.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1300<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1300.close();
      } 
                                  if(iteration_no == 1500)
      {std::ofstream alpha_field_profiles_1500("alpha_1500.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1500<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1500.close();
      } 
                                  if(iteration_no == 1700)
      {std::ofstream alpha_field_profiles_1700("alpha_1700.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_1700<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_1700.close();
      } 
                                  if(iteration_no == 2000)
      {std::ofstream alpha_field_profiles_2000("alpha_2000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_2000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_2000.close();
      } 

                                        if(iteration_no == 3000)
      {std::ofstream alpha_field_profiles_3000("alpha_3000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_3000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_3000.close();
      } 

                                        if(iteration_no == 4000)
      {std::ofstream alpha_field_profiles_4000("alpha_4000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_4000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_4000.close();
      } 

                                              if(iteration_no == 5000)
      {std::ofstream alpha_field_profiles_5000("alpha_5000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_5000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_5000.close();
      } 

                                                    if(iteration_no == 6000)
      {std::ofstream alpha_field_profiles_6000("alpha_6000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_6000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_6000.close();
      } 
                                                          if(iteration_no == 6500)
      {std::ofstream alpha_field_profiles_6500("alpha_6500.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_6500<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_6500.close();
      } 

                                                          if(iteration_no == 7000)
      {std::ofstream alpha_field_profiles_7000("alpha_7000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_7000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_7000.close();
      } 

      if(iteration_no == 7500)
      {std::ofstream alpha_field_profiles_7500("alpha_7500.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_7500<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_7500.close();
      } 

                                                          if(iteration_no == 8000)
      {std::ofstream alpha_field_profiles_8000("alpha_8000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_8000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_8000.close();
      } 

                                                          if(iteration_no == 8500)
      {std::ofstream alpha_field_profiles_8500("alpha_8500.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_8500<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_8500.close();
      } 

                                                          if(iteration_no == 9000)
      {std::ofstream alpha_field_profiles_9000("alpha_9000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_9000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_9000.close();
      } 

                                                          if(iteration_no == 9500)
      {std::ofstream alpha_field_profiles_9500("alpha_9500.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_9500<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_9500.close();
      } 

                                                                if(iteration_no == 10000)
      {std::ofstream alpha_field_profiles_10000("alpha_10000.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[j + i*x_distance.size()].get_cell_volume();
                alpha_field_profiles_10000<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<alpha.scalar_field_values(j + i*x_distance.size())<<std::endl;
              }
         }
         alpha_field_profiles_10000.close();
      } 
        
        // alpha_field_profiles_10.close();
        // alpha_field_profiles_20.close();
        // alpha_field_profiles_30.close();
        // alpha_field_profiles_40.close();
        // alpha_field_profiles_50.close();
        // alpha_field_profiles_60.close();
        // alpha_field_profiles_70.close();
        // alpha_field_profiles_80.close();
        // alpha_field_profiles_90.close();
        // alpha_field_profiles_100.close();
        // alpha_field_profiles_110.close();
        // alpha_field_profiles_120.close();
        // alpha_field_profiles_130.close();
        // alpha_field_profiles_140.close();
        // alpha_field_profiles_150.close();
        // alpha_field_profiles_160.close();
        // alpha_field_profiles_170.close();
        // alpha_field_profiles_180.close();
        // alpha_field_profiles_190.close();
        // alpha_field_profiles_200.close();
        // alpha_field_profiles_210.close();
        // alpha_field_profiles_220.close();
        // alpha_field_profiles_230.close();
        // alpha_field_profiles_240.close();  
        // alpha_field_profiles_250.close();
        // alpha_field_profiles_260.close();
        // alpha_field_profiles_270.close();
        // alpha_field_profiles_280.close();
        // alpha_field_profiles_290.close();
        // alpha_field_profiles_300.close();
} 


void PV_coupling::alpha_output_scalar_matrix_coefficients_to_file(double total_cells)
{
     std::ofstream alpha_convection_coeffs_a("alpha_convection_coeffs_a.txt");
     std::ofstream alpha_convection_coeffs_b("alpha_convection_coeffs_b.txt");
     std::ofstream alpha_roc_coeffs_a("alpha_roc_coeffs_a.txt");
     std::ofstream alpha_roc_coeffs_b("alpha_roc_coeffs_b.txt");
     std::ofstream alpha_source_coeffs_b("alpha_source_coeffs_b.txt");
     std::ofstream alpha_combined_coeffs_a("alpha_combined_coeffs_a.txt");
     std::ofstream alpha_combined_coeffs_b("alpha_combined_coeffs_b.txt");
    
     for(int i = 0; i<total_cells; i++)
     {
        for(int j=0; j<total_cells; j++)
              {
                if(j < total_cells - 1)
                {
                 alpha_convection_coeffs_a<<alpha_a_matrix_convection(i,j)<<"\t";
                 alpha_roc_coeffs_a<<alpha_a_matrix_rate_of_change(i,j)<<"\t";
                 alpha_combined_coeffs_a<<alpha_a_matrix_combined(i,j)<<"\t";

                }

                else
                {
                 alpha_convection_coeffs_a<<alpha_a_matrix_convection(i,j)<<"\n";
                 alpha_roc_coeffs_a<<alpha_a_matrix_rate_of_change(i,j)<<"\n";
                 alpha_combined_coeffs_a<<alpha_a_matrix_combined(i,j)<<"\n";

                }
              }
     }

     for(int k=0; k<total_cells; k++)
     {
                alpha_convection_coeffs_b<<alpha_b_vector_convection(k)<<"\n";
                alpha_roc_coeffs_b<<alpha_b_vector_rate_of_change(k)<<"\n";
                alpha_source_coeffs_b<<alpha_b_vector_source(k)<<"\n";
                alpha_combined_coeffs_b<<alpha_b_vector_combined(k)<<"\n";

     }

  alpha_convection_coeffs_a.close();
   alpha_convection_coeffs_b.close();
   alpha_roc_coeffs_a.close();
   alpha_roc_coeffs_b.close();
   alpha_combined_coeffs_a.close();
   alpha_combined_coeffs_b.close();
   alpha_source_coeffs_b.close();
}

