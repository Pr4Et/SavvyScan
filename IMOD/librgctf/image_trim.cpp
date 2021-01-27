#include "core_headers.h"


/* IMOD
wxMutex Image::s_mutexProtectingFFTW;*/

void Image::SetupInitialValues()
{
	logical_x_dimension = 0;
	logical_y_dimension = 0;
	logical_z_dimension = 0;

	is_in_real_space = true;
	object_is_centred_in_box = true;

	physical_upper_bound_complex_x = 0;
	physical_upper_bound_complex_y = 0;
	physical_upper_bound_complex_z = 0;

	physical_address_of_box_center_x = 0;
	physical_address_of_box_center_y = 0;
	physical_address_of_box_center_z = 0;

	//physical_index_of_first_negative_frequency_x = 0;
	physical_index_of_first_negative_frequency_y = 0;
	physical_index_of_first_negative_frequency_z = 0;

	fourier_voxel_size_x = 0.0;
	fourier_voxel_size_y = 0.0;
	fourier_voxel_size_z = 0.0;

	logical_upper_bound_complex_x = 0;
	logical_upper_bound_complex_y = 0;
	logical_upper_bound_complex_z = 0;

	logical_lower_bound_complex_x = 0;
	logical_lower_bound_complex_y = 0;
	logical_lower_bound_complex_z = 0;

	logical_upper_bound_real_x = 0;
	logical_upper_bound_real_y = 0;
	logical_upper_bound_real_z = 0;

	logical_lower_bound_real_x = 0;
	logical_lower_bound_real_y = 0;
	logical_lower_bound_real_z = 0;

	insert_into_which_reconstruction = 0;
	real_values = NULL;
	complex_values = NULL;

	is_in_memory = false;
	real_memory_allocated = 0;

	/* IMOD
  plan_fwd = NULL;
  plan_bwd = NULL;

  planned = false;*/

	padding_jump_value = 0;
	image_memory_should_not_be_deallocated = false;
}

Image::Image()
{
	SetupInitialValues();
}

Image::Image( const Image &other_image) // copy constructor
{

	SetupInitialValues();
	MyDebugPrint("Warning: copying an image object");
	*this = other_image;
}

Image::~Image()
{
	Deallocate();
}


int Image::ReturnSmallestLogicalDimension()
{
	if (logical_z_dimension == 1)
	{
		return std::min(logical_x_dimension, logical_y_dimension);
	}
	else
	{
		int temp_int;
		temp_int = std::min(logical_x_dimension, logical_y_dimension);
		return std::min(temp_int, logical_z_dimension);
	}
}


int Image::ReturnLargestLogicalDimension()
{
	if (logical_z_dimension == 1)
	{
		return std::max(logical_x_dimension, logical_y_dimension);
	}
	else
	{
		int temp_int;
		temp_int = std::max(logical_x_dimension, logical_y_dimension);
		return std::max(temp_int, logical_z_dimension);
	}
}

void Image::MultiplyPixelWise(Image &other_image)
{
	MyDebugAssertTrue(is_in_memory, "Image memory not allocated");
	MyDebugAssertTrue(other_image.is_in_memory, "Other image memory not allocated");
	MyDebugAssertTrue(is_in_real_space == other_image.is_in_real_space, "Both images need to be in same space");
	MyDebugAssertTrue(HasSameDimensionsAs(&other_image),"Images do not have same dimensions");

	int i;
	long pixel_counter;

	if (is_in_real_space)
	{
		for (pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
		{
			real_values[pixel_counter] *= other_image.real_values[pixel_counter];
		}
	}
	else
	{
		// TODO: add MKL implementation (see EulerSearch::Run for a similar example)
		for (pixel_counter = 0; pixel_counter < real_memory_allocated / 2; pixel_counter++)
		{
			complex_values[pixel_counter] *= other_image.complex_values[pixel_counter];
		}
	}
}


void Image::CircleMask(float wanted_mask_radius, bool invert)
{
	MyDebugAssertTrue(is_in_real_space,"Image not in real space");
	MyDebugAssertTrue(object_is_centred_in_box,"Object not centered in box");

	long pixel_counter;
	int i,j,k;
	float x,y,z;
	float distance_from_center_squared;
	const float wanted_mask_radius_squared = powf(wanted_mask_radius,2);
	double average_value = 0.0;
	long number_of_pixels = 0;

	pixel_counter = 0;
	for (k = 0; k < logical_z_dimension; k++)
	{
		z = powf(k - physical_address_of_box_center_z, 2);

		for (j = 0; j < logical_y_dimension; j++)
		{
			y = powf(j - physical_address_of_box_center_y, 2);

			for (i = 0; i < logical_x_dimension; i++)
			{
				x = powf(i - physical_address_of_box_center_x, 2);

				distance_from_center_squared = x + y + z;

				if (abs(distance_from_center_squared-wanted_mask_radius_squared) <= 4.0)
				{
					number_of_pixels++;
					average_value += real_values[pixel_counter];
				}

				pixel_counter++;

			}
			pixel_counter += padding_jump_value;
		}
	}

	// Now we know what value to mask with
	average_value /= float(number_of_pixels);

	// Let's mask
	pixel_counter = 0;
	for (k = 0; k < logical_z_dimension; k++)
	{
		z = powf(k - physical_address_of_box_center_z, 2);

		for (j = 0; j < logical_y_dimension; j++)
		{
			y = powf(j - physical_address_of_box_center_y, 2);

			for (i = 0; i < logical_x_dimension; i++)
			{
				x = powf(i - physical_address_of_box_center_x, 2);

				distance_from_center_squared = x + y + z;

				if (invert)
				{
					if ( distance_from_center_squared <= wanted_mask_radius_squared)
					{
						real_values[pixel_counter] = average_value;
					}
				}
				else
				{
					if ( distance_from_center_squared > wanted_mask_radius_squared)
					{
						real_values[pixel_counter] = average_value;
					}
				}

				pixel_counter++;

			}
			pixel_counter += padding_jump_value;
		}
	}


}

/* Had to import from image.cpp, it was not in a ctffind4 section */
float Image::CosineMask(float wanted_mask_radius, float wanted_mask_edge, bool invert, bool force_mask_value, float wanted_mask_value)
{
//	MyDebugAssertTrue(! is_in_real_space || object_is_centred_in_box, "Image in real space but not centered");
	MyDebugAssertTrue(wanted_mask_edge > 0, "Edge width too small");

	int i;
	int j;
	int k;
	int ii;
	int jj;
	int kk;
	long number_of_pixels;

	float x;
	float y;
	float z;

	long pixel_counter = 0;

	float distance_from_center;
	float mask_radius_plus_edge;
	float distance_from_center_squared;
	float mask_radius;
	float mask_radius_squared;
	float mask_radius_plus_edge_squared;
	float edge;
	double pixel_sum;

	float frequency;
	float frequency_squared;

	double mask_volume = 0.0;

	mask_radius = wanted_mask_radius - wanted_mask_edge / 2;
	if (mask_radius < 0.0) mask_radius = 0.0;
	mask_radius_plus_edge = mask_radius + wanted_mask_edge;

	mask_radius_squared = powf(mask_radius, 2);
	mask_radius_plus_edge_squared = powf(mask_radius_plus_edge, 2);

	pixel_sum = 0.0;
	number_of_pixels = 0;
	if (is_in_real_space && object_is_centred_in_box)
	{
		if (force_mask_value)
		{
			pixel_sum = wanted_mask_value;
		}
		else
		{
			for (k = 0; k < logical_z_dimension; k++)
			{
				z = powf(k - physical_address_of_box_center_z, 2);

				for (j = 0; j < logical_y_dimension; j++)
				{
					y = powf(j - physical_address_of_box_center_y, 2);

					for (i = 0; i < logical_x_dimension; i++)
					{
						x = powf(i - physical_address_of_box_center_x, 2);

						distance_from_center_squared = x + y + z;

						if (distance_from_center_squared >= mask_radius_squared && distance_from_center_squared <= mask_radius_plus_edge_squared)
						{
							pixel_sum += real_values[pixel_counter];
							number_of_pixels++;
						}
						pixel_counter++;
					}
					pixel_counter += padding_jump_value;
				}
			}
			pixel_sum /= number_of_pixels;
		}

		pixel_counter = 0.0;
		for (k = 0; k < logical_z_dimension; k++)
		{
			z = powf(k - physical_address_of_box_center_z, 2);

			for (j = 0; j < logical_y_dimension; j++)
			{
				y = powf(j - physical_address_of_box_center_y, 2);

				for (i = 0; i < logical_x_dimension; i++)
				{
					x = powf(i - physical_address_of_box_center_x, 2);

					distance_from_center_squared = x + y + z;

					if (distance_from_center_squared >= mask_radius_squared && distance_from_center_squared <= mask_radius_plus_edge_squared)
					{
						distance_from_center = sqrtf(distance_from_center_squared);
						edge = (1.0 + cosf(PI * (distance_from_center - mask_radius) / wanted_mask_edge)) / 2.0;
						if (invert)
						{
							real_values[pixel_counter] = real_values[pixel_counter] * (1.0 - edge) + edge * pixel_sum;
							mask_volume += powf(1.0 - edge,2);
						}
						else
						{
							real_values[pixel_counter] = real_values[pixel_counter] * edge + (1.0 - edge) * pixel_sum;
							mask_volume += powf(edge,2);
						}
					}
					else
					if (invert)
					{
						if (distance_from_center_squared <= mask_radius_squared)
						{
							real_values[pixel_counter] = pixel_sum;
						}
						else
						{
							mask_volume += 1.0;
						}
					}
					else
					{
						if (distance_from_center_squared >= mask_radius_plus_edge_squared)
						{
							real_values[pixel_counter] = pixel_sum;
						}
						else
						{
							mask_volume += 1.0;
						}
					}

					pixel_counter++;
				}
				pixel_counter += padding_jump_value;
			}
		}
	}
	else
	if (is_in_real_space)
	{
		if (force_mask_value)
		{
			pixel_sum = wanted_mask_value;
		}
		else
		{
			for (k = 0; k < logical_z_dimension; k++)
			{
				kk = k;
				if (kk >= physical_address_of_box_center_z) kk -= logical_z_dimension;
				z = powf(kk, 2);

				for (j = 0; j < logical_y_dimension; j++)
				{
					jj = j;
					if (jj >= physical_address_of_box_center_y) jj -= logical_y_dimension;
					y = powf(jj, 2);

					for (i = 0; i < logical_x_dimension; i++)
					{
						ii = i;
						if (ii >= physical_address_of_box_center_x) ii -= logical_x_dimension;
						x = powf(ii, 2);

						distance_from_center_squared = x + y + z;

						if (distance_from_center_squared >= mask_radius_squared && distance_from_center_squared <= mask_radius_plus_edge_squared)
						{
							pixel_sum += real_values[pixel_counter];
							number_of_pixels++;
						}
						pixel_counter++;
					}
					pixel_counter += padding_jump_value;
				}
			}
			pixel_sum /= number_of_pixels;
		}

		pixel_counter = 0.0;
		for (k = 0; k < logical_z_dimension; k++)
		{
			kk = k;
			if (kk >= physical_address_of_box_center_z) kk -= logical_z_dimension;
			z = powf(kk, 2);

			for (j = 0; j < logical_y_dimension; j++)
			{
				jj = j;
				if (jj >= physical_address_of_box_center_y) jj -= logical_y_dimension;
				y = powf(jj, 2);

				for (i = 0; i < logical_x_dimension; i++)
				{
					ii = i;
					if (ii >= physical_address_of_box_center_x) ii -= logical_x_dimension;
					x = powf(ii, 2);

					distance_from_center_squared = x + y + z;

					if (distance_from_center_squared >= mask_radius_squared && distance_from_center_squared <= mask_radius_plus_edge_squared)
					{
						distance_from_center = sqrtf(distance_from_center_squared);
						edge = (1.0 + cosf(PI * (distance_from_center - mask_radius) / wanted_mask_edge)) / 2.0;
						real_values[pixel_counter] = real_values[pixel_counter] * edge + (1.0 - edge) * pixel_sum;
						mask_volume += powf(edge,2);
					}
					else
						if (distance_from_center_squared >= mask_radius_plus_edge_squared) real_values[pixel_counter] = pixel_sum;
					else
					{
						mask_volume += 1.0;
					}

					pixel_counter++;
				}
				pixel_counter += padding_jump_value;
			}
		}
	}
	else
	{
		for (k = 0; k <= physical_upper_bound_complex_z; k++)
		{
			z = powf(ReturnFourierLogicalCoordGivenPhysicalCoord_Z(k) * fourier_voxel_size_z, 2);

			for (j = 0; j <= physical_upper_bound_complex_y; j++)
			{
				y = powf(ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * fourier_voxel_size_y, 2);

				for (i = 0; i <= physical_upper_bound_complex_x; i++)
				{
					x = powf(i * fourier_voxel_size_x, 2);

					// compute squared radius, in units of reciprocal pixels

					frequency_squared = x + y + z;

					if (frequency_squared >= mask_radius_squared && frequency_squared <= mask_radius_plus_edge_squared)
					{
						frequency = sqrtf(frequency_squared);
						edge = (1.0 + cosf(PI * (frequency - mask_radius) / wanted_mask_edge)) / 2.0;
						if (invert)
						{
							complex_values[pixel_counter] *= (1.0 - edge);
						}
						else
						{
							complex_values[pixel_counter] *= edge;
						}
					}
					if (invert)
					{
						if (frequency_squared <= mask_radius_squared) complex_values[pixel_counter] = 0.0f + I * 0.0f;
					}
					else
					{
						if (frequency_squared >= mask_radius_plus_edge_squared) complex_values[pixel_counter] = 0.0f + I * 0.0f;
					}

					pixel_counter++;
				}
			}
		}
	}
	
	return float(mask_volume);
}

Image & Image::operator = (const Image &other_image)
{
	*this = &other_image;
	return *this;
}


Image & Image::operator = (const Image *other_image)
{
   // Check for self assignment
   if(this != other_image)
   {
		MyDebugAssertTrue(other_image->is_in_memory, "Other image Memory not allocated");

		if (is_in_memory == true)
		{

			if (logical_x_dimension != other_image->logical_x_dimension || logical_y_dimension != other_image->logical_y_dimension || logical_z_dimension != other_image->logical_z_dimension)
			{
				Deallocate();
				Allocate(other_image->logical_x_dimension, other_image->logical_y_dimension, other_image->logical_z_dimension, other_image->is_in_real_space);
			}
		}
		else
		{
			Allocate(other_image->logical_x_dimension, other_image->logical_y_dimension, other_image->logical_z_dimension, other_image->is_in_real_space);
		}

		// by here the memory allocation should be ok..

		is_in_real_space = other_image->is_in_real_space;
		object_is_centred_in_box = other_image->object_is_centred_in_box;

		for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
		{
			real_values[pixel_counter] = other_image->real_values[pixel_counter];
		}
   }

   return *this;
}



//!>  \brief  Deallocate all the memory.  The idea is to keep this safe in the case that something isn't
//    allocated, so it can always be safely called.  I.e. it should check whatever it is deallocating is
//    in fact allocated.
//


void Image::Deallocate()
{
	if (is_in_memory == true && image_memory_should_not_be_deallocated == false)
	{
		/* IMOD
    fftwf_free(real_values);*/
		free(real_values);
		is_in_memory = false;
	}

  /* IMOD
	if (planned == true)
	{
		wxMutexLocker lock(s_mutexProtectingFFTW); // the mutex will be unlocked when this object is destroyed (when it goes out of scope)
    	MyDebugAssertTrue(lock.IsOk(),"Mute locking failed");
		fftwf_destroy_plan(plan_fwd);
		fftwf_destroy_plan(plan_bwd);
		planned = false;
    }*/

}

//!>  \brief  Allocate memory for the Image object.
//
//  If the object is already allocated with correct dimensions, nothing happens. Otherwise, object is deallocated first.

void Image::Allocate(int wanted_x_size, int wanted_y_size, int wanted_z_size, bool should_be_in_real_space)
{

	MyDebugAssertTrue(wanted_x_size > 0 && wanted_y_size > 0 && wanted_z_size > 0,"Bad dimensions: %i %i %i\n",wanted_x_size,wanted_y_size,wanted_z_size);

	// check to see if we need to do anything?

	if (is_in_memory == true)
	{
		is_in_real_space = should_be_in_real_space;

		if (wanted_x_size == logical_x_dimension && wanted_y_size == logical_y_dimension && wanted_z_size == logical_z_dimension)
		{
			// everything is already done..

			is_in_real_space = should_be_in_real_space;
//			wxPrintf("returning\n");

			return;
		}
		else
		Deallocate();
	}

	// if we got here we need to do allocation..

	SetLogicalDimensions(wanted_x_size, wanted_y_size, wanted_z_size);
	is_in_real_space = should_be_in_real_space;

	// first_x_dimension
	if (IsEven(wanted_x_size) == true) real_memory_allocated =  wanted_x_size / 2 + 1;
	else real_memory_allocated = (wanted_x_size - 1) / 2 + 1;

	real_memory_allocated *= wanted_y_size * wanted_z_size; // other dimensions
	real_memory_allocated *= 2; // room for complex

	/* IMOD  switch to regular malloc
  real_values = (float *) fftwf_malloc(sizeof(float) * real_memory_allocated);*/
	real_values = (float *) malloc(sizeof(float) * real_memory_allocated);
	complex_values = (std::complex<float>*) real_values;  // Set the complex_values to point at the newly allocated real values;

	is_in_memory = true;

	// Update addresses etc..

    UpdateLoopingAndAddressing();

    // Prepare the plans for FFTW

    /* IMOD
    if (planned == false)
    {
    	wxMutexLocker lock(s_mutexProtectingFFTW); // the mutex will be unlocked when this object is destroyed (when it goes out of scope)
    	MyDebugAssertTrue(lock.IsOk(),"Mute locking failed");
    	if (logical_z_dimension > 1)
    	{
    		plan_fwd = fftwf_plan_dft_r2c_3d(logical_z_dimension, logical_y_dimension, logical_x_dimension, real_values, reinterpret_cast<fftwf_complex*>(complex_values), FFTW_ESTIMATE);
    		plan_bwd = fftwf_plan_dft_c2r_3d(logical_z_dimension, logical_y_dimension, logical_x_dimension, reinterpret_cast<fftwf_complex*>(complex_values), real_values, FFTW_ESTIMATE);
    	}
    	else
    	{
    		plan_fwd = fftwf_plan_dft_r2c_2d(logical_y_dimension, logical_x_dimension, real_values, reinterpret_cast<fftwf_complex*>(complex_values), FFTW_ESTIMATE);
    	    plan_bwd = fftwf_plan_dft_c2r_2d(logical_y_dimension, logical_x_dimension, reinterpret_cast<fftwf_complex*>(complex_values), real_values, FFTW_ESTIMATE);

    	}

    	planned = true;
    }*/

    // set the loop junk value..

	if (IsEven(logical_x_dimension) == true) padding_jump_value = 2;
	else padding_jump_value = 1;

	//

	number_of_real_space_pixels = long(logical_x_dimension) * long(logical_y_dimension) * long(logical_z_dimension);
	ft_normalization_factor = 1.0 / sqrtf(float(number_of_real_space_pixels));
}

//!>  \brief  Allocate memory for the Image object.
//
//  Overloaded version of allocate to cover the supplying just 2 dimensions along with the should_be_in_real_space bool.

void Image::Allocate(int wanted_x_size, int wanted_y_size, bool should_be_in_real_space)
{
	Allocate(wanted_x_size, wanted_y_size, 1, should_be_in_real_space);
}


void Image::AllocateAsPointingToSliceIn3D(Image *wanted3d, long wanted_slice)
{
	Deallocate();
	is_in_real_space = wanted3d->is_in_real_space;

	// if we got here we need to do allocation..

	SetLogicalDimensions(wanted3d->logical_x_dimension, wanted3d->logical_y_dimension, 1);

	// we are not actually allocating, we are pointing..

	long bytes_in_slice = wanted3d->real_memory_allocated / wanted3d->logical_z_dimension;

	image_memory_should_not_be_deallocated = true;
	is_in_memory = true; // kind of a lie
	real_memory_allocated = bytes_in_slice; // kind of a lie

	real_values = wanted3d->real_values + (bytes_in_slice * (wanted_slice - 1)); // point to the 3d..
	complex_values = (std::complex<float>*) real_values;  // Set the complex_values to point at the newly allocated real values;

	// Update addresses etc..

    UpdateLoopingAndAddressing();

    // Prepare the plans for FFTW

    /* IMOD
    if (planned == false)
    {
    	wxMutexLocker lock(s_mutexProtectingFFTW); // the mutex will be unlocked when this object is destroyed (when it goes out of scope)
        MyDebugAssertTrue(lock.IsOk(),"Mute locking failed");

    	if (logical_z_dimension > 1)
    	{
    		plan_fwd = fftwf_plan_dft_r2c_3d(logical_z_dimension, logical_y_dimension, logical_x_dimension, real_values, reinterpret_cast<fftwf_complex*>(complex_values), FFTW_ESTIMATE);
    		plan_bwd = fftwf_plan_dft_c2r_3d(logical_z_dimension, logical_y_dimension, logical_x_dimension, reinterpret_cast<fftwf_complex*>(complex_values), real_values, FFTW_ESTIMATE);
    	}
    	else
    	{
    		plan_fwd = fftwf_plan_dft_r2c_2d(logical_y_dimension, logical_x_dimension, real_values, reinterpret_cast<fftwf_complex*>(complex_values), FFTW_ESTIMATE);
    	    plan_bwd = fftwf_plan_dft_c2r_2d(logical_y_dimension, logical_x_dimension, reinterpret_cast<fftwf_complex*>(complex_values), real_values, FFTW_ESTIMATE);

    	}

    	planned = true;
    }*/

    // set the loop jump value..

	if (IsEven(logical_x_dimension) == true) padding_jump_value = 2;
	else padding_jump_value = 1;

	//
	number_of_real_space_pixels = long(logical_x_dimension) * long(logical_y_dimension) * long(logical_z_dimension);
	ft_normalization_factor = 1.0 / sqrtf(float(number_of_real_space_pixels));

}

//!>  \brief  Change the logical dimensions of an image and update all related values

void Image::SetLogicalDimensions(int wanted_x_size, int wanted_y_size, int wanted_z_size)
{
	logical_x_dimension = wanted_x_size;
	logical_y_dimension = wanted_y_size;
	logical_z_dimension = wanted_z_size;
}

//!>  \brief  Update all properties related to looping & addressing in real & Fourier space, given the current logical dimensions.

void Image::UpdateLoopingAndAddressing()
{

	physical_upper_bound_complex_x = logical_x_dimension / 2;
	physical_upper_bound_complex_y = logical_y_dimension - 1;
	physical_upper_bound_complex_z = logical_z_dimension - 1;

	UpdatePhysicalAddressOfBoxCenter();

	//physical_index_of_first_negative_frequency_x = logical_x_dimension / 2 + 1;
	if (IsEven(logical_y_dimension) == true)
	{
		physical_index_of_first_negative_frequency_y = logical_y_dimension / 2;
	}
	else
	{
		physical_index_of_first_negative_frequency_y = logical_y_dimension / 2 + 1;
	}

	if (IsEven(logical_z_dimension) == true)
	{
		physical_index_of_first_negative_frequency_z = logical_z_dimension / 2;
	}
	else
	{
		physical_index_of_first_negative_frequency_z = logical_z_dimension / 2 + 1;
	}


    // Update the Fourier voxel size

	fourier_voxel_size_x = 1.0 / double(logical_x_dimension);
	fourier_voxel_size_y = 1.0 / double(logical_y_dimension);
	fourier_voxel_size_z = 1.0 / double(logical_z_dimension);

	// Logical bounds
	if (IsEven(logical_x_dimension) == true)
	{
		logical_lower_bound_complex_x = -logical_x_dimension / 2;
		logical_upper_bound_complex_x =  logical_x_dimension / 2;
	    logical_lower_bound_real_x    = -logical_x_dimension / 2;
	    logical_upper_bound_real_x    =  logical_x_dimension / 2 - 1;
	}
	else
	{
		logical_lower_bound_complex_x = -(logical_x_dimension-1) / 2;
		logical_upper_bound_complex_x =  (logical_x_dimension-1) / 2;
		logical_lower_bound_real_x    = -(logical_x_dimension-1) / 2;
		logical_upper_bound_real_x    =  (logical_x_dimension-1) / 2;
	}


	if (IsEven(logical_y_dimension) == true)
	{
	    logical_lower_bound_complex_y = -logical_y_dimension / 2;
	    logical_upper_bound_complex_y =  logical_y_dimension / 2 - 1;
	    logical_lower_bound_real_y    = -logical_y_dimension / 2;
	    logical_upper_bound_real_y    =  logical_y_dimension / 2 - 1;
	}
	else
	{
	    logical_lower_bound_complex_y = -(logical_y_dimension-1) / 2;
	    logical_upper_bound_complex_y =  (logical_y_dimension-1) / 2;
	    logical_lower_bound_real_y    = -(logical_y_dimension-1) / 2;
	    logical_upper_bound_real_y    =  (logical_y_dimension-1) / 2;
	}

	if (IsEven(logical_z_dimension) == true)
	{
		logical_lower_bound_complex_z = -logical_z_dimension / 2;
		logical_upper_bound_complex_z =  logical_z_dimension / 2 - 1;
		logical_lower_bound_real_z    = -logical_z_dimension / 2;
		logical_upper_bound_real_z    =  logical_z_dimension / 2 - 1;

	}
	else
	{
		logical_lower_bound_complex_z = -(logical_z_dimension - 1) / 2;
		logical_upper_bound_complex_z =  (logical_z_dimension - 1) / 2;
		logical_lower_bound_real_z    = -(logical_z_dimension - 1) / 2;
		logical_upper_bound_real_z    =  (logical_z_dimension - 1) / 2;
	}
}

//!>  \brief  Returns the physical address of the image origin

void Image::UpdatePhysicalAddressOfBoxCenter()
{
	/*
    if (IsEven(logical_x_dimension)) physical_address_of_box_center_x = logical_x_dimension / 2;
    else physical_address_of_box_center_x = (logical_x_dimension - 1) / 2;

    if (IsEven(logical_y_dimension)) physical_address_of_box_center_y = logical_y_dimension / 2;
    else physical_address_of_box_center_y = (logical_y_dimension - 1) / 2;

    if (IsEven(logical_z_dimension)) physical_address_of_box_center_z = logical_z_dimension / 2;
    else physical_address_of_box_center_z = (logical_z_dimension - 1) / 2;
*/
	physical_address_of_box_center_x = logical_x_dimension / 2;
	physical_address_of_box_center_y = logical_y_dimension / 2;
	physical_address_of_box_center_z = logical_z_dimension / 2;
}


// Work out whether a given Fourier component has a Hermitian mate which is also described explicitely by the FFTW
bool Image::FourierComponentHasExplicitHermitianMate(int physical_index_x, int physical_index_y, int physical_index_z)
{
	bool explicit_mate;

	explicit_mate = physical_index_x == 0 && ! ( physical_index_y == 0 && physical_index_z == 0);

	// We assume that the Y dimension is the non-flat one
	if (IsEven(logical_y_dimension))
	{
		explicit_mate = explicit_mate && physical_index_y != physical_index_of_first_negative_frequency_y-1;
	}

	if (logical_z_dimension > 1)
	{
		if (IsEven(logical_z_dimension))
		{
			explicit_mate = explicit_mate && physical_index_z != physical_index_of_first_negative_frequency_z - 1;
		}
	}

	return explicit_mate;
}

// Work out whether a given Fourier component is a (redundant) Hermitian mate which is described explicitely by the FFTW but
// shouldn't be counted in statistics as an independent Fourier component
bool Image::FourierComponentIsExplicitHermitianMate(int physical_index_x, int physical_index_y, int physical_index_z)
{
	bool explicit_mate = physical_index_x == 0 && (physical_index_y >= physical_index_of_first_negative_frequency_y || physical_index_z >= physical_index_of_first_negative_frequency_z);

	return explicit_mate;
}


//!> \brief   Apply a forward FT to the Image object. The FT is scaled.
//   The DC component is at (self%DIM(1)/2+1,self%DIM(2)/2+1,self%DIM(3)/2+1) (assuming even dimensions) or at (1,1,1) by default.
//
//
//   For details on FFTW, see http://www.fftw.org/
//   A helpful page for understanding the output format: http://www.dsprelated.com/showmessage/102675/1.php
//   A helpful page to learn about vectorization and FFTW benchmarking: http://www.dsprelated.com/showmessage/76779/1.php
//   \todo   Check this: http://objectmix.com/fortran/371439-ifort-openmp-fftw-problem.html
//
//
//
// 	A note on scaling: by default, we divide by N, the number of pixels. This ensures that after we do an inverse FT (without further scaling),
//	we will return to our original values. However, it means that while in Fourier space, the amplitudes are too high, by a factor of sqrt(N),
//  such that, for example, Parserval's theorem is not satisfied.
void Image::ForwardFFT(bool should_scale)
{

	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Image already in Fourier space");
	/* IMOD
  MyDebugAssertTrue(planned, "FFT's not planned");

	fftwf_execute_dft_r2c(plan_fwd, real_values, reinterpret_cast<fftwf_complex*>(complex_values));

	if (should_scale)
	{
		DivideByConstant(float(number_of_real_space_pixels));
	}

	// Set the image type

	is_in_real_space = false;*/
}

void Image::BackwardFFT()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertFalse(is_in_real_space, "Image already in real space");

	/* IMOD
  fftwf_execute_dft_c2r(plan_bwd, reinterpret_cast<fftwf_complex*>(complex_values), real_values);

    // Set the image type

    is_in_real_space = true; */
}

//!> \brief Divide all voxels by a constant value (this is actually done as a multiplication by the inverse)

void Image::DivideByConstant(float constant_to_divide_by)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	float inverse = 1. / constant_to_divide_by;
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] *= inverse;
	}
}

void Image::MultiplyAddConstant(float constant_to_multiply_by, float constant_to_add)
{
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] = real_values[pixel_counter] * constant_to_multiply_by + constant_to_add;
	}
}

void Image::AddMultiplyConstant(float constant_to_add, float constant_to_multiply_by)
{
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] = (real_values[pixel_counter] + constant_to_add) * constant_to_multiply_by;
	}
}

void Image::AddMultiplyAddConstant(float first_constant_to_add, float constant_to_multiply_by, float second_constant_to_add)
{
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] = (real_values[pixel_counter] + first_constant_to_add) * constant_to_multiply_by + second_constant_to_add;
	}
}

//!> \brief Multiply all voxels by a constant value

//inline
void Image::MultiplyByConstant(float constant_to_multiply_by)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] *= constant_to_multiply_by;
	}
}

void Image::TakeReciprocalRealValues(float zeros_become)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		if (real_values[pixel_counter] != 0.0) real_values[pixel_counter] = 1.0 / real_values[pixel_counter];
		else real_values[pixel_counter] = zeros_become;
	}

}

void Image::InvertRealValues()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] = - real_values[pixel_counter];
	}
}

void Image::SquareRealValues()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Must be in real space to square real values");
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter ++ )
	{
		real_values[pixel_counter] *= real_values[pixel_counter];
	}
}

void Image::ExponentiateRealValues()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Must be in real space to square real values");
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter ++ )
	{
		real_values[pixel_counter] = exp(real_values[pixel_counter]);
	}
}

void Image::SquareRootRealValues()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Must be in real space to square real values");
	MyDebugAssertFalse(HasNegativeRealValue(),"Image has negative value(s). Cannot compute square root.\n");
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter ++ )
	{
		real_values[pixel_counter] = sqrtf(real_values[pixel_counter]);
	}
}

bool Image::IsConstant()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		if (real_values[pixel_counter] != real_values[0]) return false;
	}
	return true;
}

bool Image::IsBinary()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Only makes sense for images in real space")

	long pixel_counter = 0;

	for ( int k = 0; k < logical_z_dimension; k ++ )
	{
		for ( int j = 0; j < logical_y_dimension; j ++ )
		{
			for ( int i = 0; i < logical_x_dimension; i ++ )
			{
				if (real_values[pixel_counter] != 0 && real_values[pixel_counter] != 1) return false;
				pixel_counter++;
			}
			pixel_counter += padding_jump_value;
		}
	}

	return true;
}

bool Image::HasNan()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	long pixel_counter = 0;
	for ( int k = 0; k < logical_z_dimension; k ++ )
	{
		for ( int j = 0; j < logical_y_dimension; j ++ )
		{
			for ( int i = 0; i < logical_x_dimension; i ++ )
			{
        /* IMOD icl 11 and maybe icc 11: use isnan instead of std::isnan */
#ifdef _WIN32
        #define isnan _isnan
#endif
				if (isnan(real_values[pixel_counter])) return true; // isnan() is part of C++11
				//if (std::isnan(real_values[pixel_counter])) return true; // std::isnan() is part of C++14
				pixel_counter ++;
			}
			pixel_counter += padding_jump_value;
		}


	}
	return false;
}

bool Image::HasNegativeRealValue()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	long pixel_counter = 0;
	for ( int k = 0; k < logical_z_dimension; k ++ )
	{
		for ( int j = 0; j < logical_y_dimension; j ++ )
		{
			for ( int i = 0; i < logical_x_dimension; i ++ )
			{
				if (real_values[pixel_counter] < 0.0) return true;
				pixel_counter ++;
			}
			pixel_counter += padding_jump_value;
		}
	}
	return false;
}

/* IMOD
void Image::ReadSlices(ImageFile *input_file, long start_slice, long end_slice)
{
	MyDebugAssertTrue(start_slice <= end_slice, "Start slice larger than end slice!");
	MyDebugAssertTrue(start_slice > 0, "Start slice is less than 0, the first slice is 1!");
	MyDebugAssertTrue(end_slice <= input_file->ReturnNumberOfSlices(), "End slice is greater than number of slices in the file!");
	MyDebugAssertTrue(input_file->IsOpen(), "Image file is not open!");


	// check the allocations..

	int number_of_slices = (end_slice - start_slice) + 1;

	if (logical_x_dimension != input_file->ReturnXSize() || logical_y_dimension != input_file->ReturnYSize() || logical_z_dimension != number_of_slices || is_in_memory == false)
	{
		Deallocate();
		Allocate(input_file->ReturnXSize(), input_file->ReturnYSize(), number_of_slices);

	}

	// We should be set up - so read in the correct values.
	// AT THE MOMENT, WE CAN ONLY READ REAL SPACE IMAGES, SO OVERRIDE THIS!!!!!

	is_in_real_space = true;
	object_is_centred_in_box = true;

	input_file->ReadSlicesFromDisk(start_slice, end_slice, real_values);

	// we need to respace this to take into account the FFTW padding..

	AddFFTWPadding();
}


//!> \brief Read a set of slices from disk (FFTW padding is done automatically)

void Image::ReadSlices(MRCFile *input_file, long start_slice, long end_slice)
{

	MyDebugAssertTrue(start_slice <= end_slice, "Start slice larger than end slice!");
	MyDebugAssertTrue(start_slice > 0, "Start slice is less than 0, the first slice is 1!");
	MyDebugAssertTrue(end_slice <= input_file->ReturnNumberOfSlices(), "End slice is greater than number of slices in the file!");
	MyDebugAssertTrue(input_file->my_file.is_open(), "MRCFile not open!");


	// check the allocations..

	int number_of_slices = (end_slice - start_slice) + 1;

	if (logical_x_dimension != input_file->ReturnXSize() || logical_y_dimension != input_file->ReturnYSize() || logical_z_dimension != number_of_slices || is_in_memory == false)
	{
		Deallocate();
		Allocate(input_file->ReturnXSize(), input_file->ReturnYSize(), number_of_slices);

	}

	// We should be set up - so read in the correct values.
	// AT THE MOMENT, WE CAN ONLY READ REAL SPACE IMAGES, SO OVERRIDE THIS!!!!!

	is_in_real_space = true;
	object_is_centred_in_box = true;

	input_file->ReadSlicesFromDisk(start_slice, end_slice, real_values);

	// we need to respace this to take into account the FFTW padding..

	AddFFTWPadding();

}

//!> \brief Read a set of slices from disk (FFTW padding is done automatically)

void Image::ReadSlices(DMFile *input_file, long start_slice, long end_slice)
{

	MyDebugAssertTrue(start_slice <= end_slice, "Start slice larger than end slice!");
	MyDebugAssertTrue(start_slice > 0, "Start slice is less than 0, the first slice is 1!");
	MyDebugAssertTrue(end_slice <= input_file->ReturnNumberOfSlices(), "End slice is greater than number of slices in the file!");
	MyDebugAssertTrue(start_slice == end_slice, "Can only read one slice at a time from DM files. Sorry.")


	// check the allocations..

	int number_of_slices = (end_slice - start_slice) + 1;

	if (logical_x_dimension != input_file->ReturnXSize() || logical_y_dimension != input_file->ReturnYSize() || logical_z_dimension != number_of_slices || is_in_memory == false)
	{
		Deallocate();
		Allocate(input_file->ReturnXSize(), input_file->ReturnYSize(), number_of_slices);

	}

	// We should be set up - so read in the correct values.
	// AT THE MOMENT, WE CAN ONLY READ REAL SPACE IMAGES, SO OVERRIDE THIS!!!!!

	is_in_real_space = true;
	object_is_centred_in_box = true;

	input_file->ReadSliceFromDisk(start_slice - 1, real_values); // DM indexes slices starting at 0

	// we need to respace this to take into account the FFTW padding..

	AddFFTWPadding();

}

//!> \brief Write a set of slices from disk (FFTW padding is done automatically)

void Image::WriteSlices(MRCFile *input_file, long start_slice, long end_slice)
{
	MyDebugAssertTrue(start_slice <= end_slice, "Start slice larger than end slice!");
	MyDebugAssertTrue(input_file->my_file.is_open(), "MRCFile not open!");

	// THIS PROBABLY NEEDS ATTENTION..

	if (start_slice == 1) // if the start slice is one, we set the header to match the image
	{
		input_file->my_header.SetDimensionsImage(logical_x_dimension,logical_y_dimension);

		if (end_slice > input_file->ReturnNumberOfSlices())
		{
			input_file->my_header.SetNumberOfImages(end_slice);
		}

		//input_file->WriteHeader();
		input_file->rewrite_header_on_close = true;
	}
	else // if the last slice is bigger than the current max number of slices, increase the max number of slices
	{
		if (end_slice > input_file->ReturnNumberOfSlices())
		{
			input_file->my_header.SetNumberOfImages(end_slice);
		}

		input_file->rewrite_header_on_close = true;
	}

	MyDebugAssertTrue(logical_x_dimension == input_file->ReturnXSize() || logical_y_dimension == input_file->ReturnYSize(), "Image dimensions (%i, %i) and file dimensions (%i, %i) differ!", logical_x_dimension, logical_y_dimension, input_file->ReturnXSize(), input_file->ReturnYSize());

	// if the image is complex.. make a temp image and transform it..

	int number_of_slices = (end_slice - start_slice) + 1;

	if (is_in_real_space == false)
	{
		Image temp_image;
		temp_image.CopyFrom(this);
		temp_image.BackwardFFT();
		temp_image.RemoveFFTWPadding();
		input_file->WriteSlicesToDisk(start_slice, end_slice, temp_image.real_values);

	}
	else // real space
	{
		RemoveFFTWPadding();
		input_file->WriteSlicesToDisk(start_slice, end_slice, real_values);
		AddFFTWPadding(); // to go back
	}
}

void Image::QuickAndDirtyWriteSlices(std::string filename, long first_slice_to_write, long last_slice_to_write)
{
	MyDebugAssertTrue(first_slice_to_write >0, "Slice is less than 1, first slice is 1");
	MRCFile output_file(filename, false);
	WriteSlices(&output_file,first_slice_to_write,last_slice_to_write);
}


void Image::QuickAndDirtyWriteSlice(std::string filename, long slice_to_write)
{
	MyDebugAssertTrue(slice_to_write >0, "Slice is less than 1, first slice is 1");
	MRCFile output_file(filename, false);
	WriteSlice(&output_file, slice_to_write);
}

void Image::QuickAndDirtyReadSlice(std::string filename, long slice_to_read)
{
	ImageFile image_file(filename);
	ReadSlice(&image_file,slice_to_read);
}
*/

/* IMOD version of write slice */
static WriteSliceType sSliceWriteFunc = NULL;

void internalSetWriteSliceFunc(WriteSliceType func)
{
  sSliceWriteFunc = func;
}

void Image::QuickAndDirtyWriteSlice(std::string filename, long slice_to_write)
{
	MyDebugAssertTrue(slice_to_write == 1, "Slice is not 1, can write only 1 slice");
  float *fdata = real_values;
  int iy;
  if (!sSliceWriteFunc)
    return;
  if (padding_jump_value > 0) {
    fdata = (float *)malloc(logical_x_dimension * logical_y_dimension * sizeof(float));
    if (!fdata)
      return;
    for (iy = 0; iy < logical_y_dimension; iy++)
      memcpy(fdata + iy * logical_x_dimension, real_values + iy * (logical_x_dimension + padding_jump_value), logical_x_dimension * sizeof(float));
  }
  sSliceWriteFunc(filename.c_str(), fdata, logical_x_dimension, logical_y_dimension);
  if (padding_jump_value > 0)
    free(fdata);
}


//!> \brief Take a contiguous set of values, and add the FFTW padding.

void Image::AddFFTWPadding()
{
	MyDebugAssertTrue(is_in_memory, "Image not allocated!");

	int x,y,z;

	long current_write_position = real_memory_allocated - (1 + padding_jump_value);
	long current_read_position = (long(logical_x_dimension) * long(logical_y_dimension) * long(logical_z_dimension)) - 1;

	for (z = 0; z < logical_z_dimension; z++)
	{
		for (y = 0; y < logical_y_dimension; y++)
		{
			for (x = 0; x < logical_x_dimension; x++)
			{
				real_values[current_write_position] = real_values[current_read_position];
				current_write_position--;
				current_read_position--;
			}

			current_write_position -= padding_jump_value;
		}
	}
}

//!> \brief Take a set of FFTW padded values, and remove the padding.

void Image::RemoveFFTWPadding()
{
	MyDebugAssertTrue(is_in_memory, "Image not allocated!");

	int x,y,z;

	long current_write_position = 0;
	long current_read_position = 0;

	for (z = 0; z < logical_z_dimension; z++)
	{
		for (y = 0; y < logical_y_dimension; y++)
		{
			for (x = 0; x < logical_x_dimension; x++)
			{
				real_values[current_write_position] = real_values[current_read_position];
				current_write_position++;
				current_read_position++;
			}

			current_read_position +=padding_jump_value;
		}
	}
}

void Image::SetToConstant(float wanted_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] = wanted_value;
	}
}

void Image::AddConstant(float wanted_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] += wanted_value;
	}
}

// DNM: Switch to looping over actual image to make valgrind happy
void Image::SetMaximumValue(float new_maximum_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i, j, k;
	long address = 0;

	for (k = 0; k < logical_z_dimension; k++)
		{
			for (j = 0; j < logical_y_dimension; j++)
			{
				for (i = 0; i < logical_x_dimension; i++)
				{
          real_values[address] = std::min(new_maximum_value,real_values[address]);
					address++;
				}
				address += padding_jump_value;
			}
		}
}

// DNM: Switch to looping over actual image to make valgrind happy
void Image::SetMinimumValue(float new_minimum_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i, j, k;
	long address = 0;

	for (k = 0; k < logical_z_dimension; k++)
		{
			for (j = 0; j < logical_y_dimension; j++)
			{
				for (i = 0; i < logical_x_dimension; i++)
				{
          real_values[address] = std::max(new_minimum_value,real_values[address]);
					address++;
				}
				address += padding_jump_value;
			}
		}
}

void Image::Binarise(float threshold_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	for (long address = 0; address < real_memory_allocated; address++)
	{
		if (real_values[address] >= threshold_value) real_values[address] = 1.0f;
		else real_values[address] = 0.0f;
	}
}


// DNM: Switch to looping over actual image to make valgrind happy
void Image::SetMinimumAndMaximumValues(float new_minimum_value, float new_maximum_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i, j, k;
	long address = 0;

	for (k = 0; k < logical_z_dimension; k++)
		{
			for (j = 0; j < logical_y_dimension; j++)
			{
				for (i = 0; i < logical_x_dimension; i++)
				{
          real_values[address] = std::max(std::min(new_maximum_value,real_values[address]),new_minimum_value);
					address++;
				}
				address += padding_jump_value;
			}
		}
}

float Image::ReturnMaximumDiagonalRadius()
{
	if (is_in_real_space)
	{
    /* IMOD icl 11 add (double) */
		return sqrt(pow((double)physical_address_of_box_center_x,2)+pow((double)physical_address_of_box_center_y,2)+pow((double)physical_address_of_box_center_z,2));
	}
	else
	{
		return sqrt(pow(logical_lower_bound_complex_x * fourier_voxel_size_x , 2) + pow(logical_lower_bound_complex_y * fourier_voxel_size_y , 2) + pow(logical_lower_bound_complex_z * fourier_voxel_size_z , 2) );
	}
}

void Image::GetMinMax(float &min_value, float &max_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Only real space supported");

	min_value = FLT_MAX;
	max_value = -FLT_MAX;

	int i, j, k;
	long address = 0;

	for (k = 0; k < logical_z_dimension; k++)
		{
			for (j = 0; j < logical_y_dimension; j++)
			{
				for (i = 0; i < logical_x_dimension; i++)
				{
					if (real_values[address] < min_value) min_value = real_values[address];
					if (real_values[address] > max_value) max_value = real_values[address];

					address++;
				}
				address += padding_jump_value;
			}
		}
}


float Image::ReturnAverageOfRealValuesOnEdges()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	double sum;
	long number_of_pixels;
	int pixel_counter;
	int line_counter;
	int plane_counter;
	long address;

	sum = 0.0;
	number_of_pixels = 0;
	address = 0;

	if (logical_z_dimension == 1)
	{
		// Two-dimensional image

		// First line
		for (pixel_counter=0; pixel_counter < logical_x_dimension; pixel_counter++)
		{
			sum += real_values[address];
			address++;
		}
		number_of_pixels += logical_x_dimension;
		address += padding_jump_value;
		// Other lines
		for (line_counter=1; line_counter < logical_y_dimension-1; line_counter++)
		{
			sum += real_values[address];
			address += logical_x_dimension-1;
			sum += real_values[address];
			address += padding_jump_value + 1;
			number_of_pixels += 2;
		}
		// Last line
		for (pixel_counter=0; pixel_counter < logical_x_dimension; pixel_counter++)
		{
			sum += real_values[address];
			address++;
		}
		number_of_pixels += logical_x_dimension;
	}
	else
	{
		// Three-dimensional volume

		// First plane
		for (line_counter=0; line_counter < logical_y_dimension; line_counter++)
		{
			for (pixel_counter=0; pixel_counter < logical_x_dimension; pixel_counter++)
			{
				sum += real_values[address];
				address++;
			}
			address += padding_jump_value;
		}
		number_of_pixels += logical_x_dimension * logical_y_dimension;
		// Other planes
		for (plane_counter = 1; plane_counter < logical_z_dimension - 1; plane_counter++)
		{
			for (line_counter=0; line_counter< logical_y_dimension; line_counter++)
			{
				if (line_counter == 0 || line_counter == logical_y_dimension-1)
				{
					// First and last line of that section
					for (pixel_counter=0; pixel_counter < logical_x_dimension; pixel_counter++)
					{
						sum += real_values[address];
						address++;
					}
					address += padding_jump_value;
					number_of_pixels += logical_x_dimension;
				}
				else
				{
					// All other lines (only count first and last pixel)
					sum += real_values[address];
					address += logical_x_dimension-1;
					sum += real_values[address];
					address += padding_jump_value + 1;
					number_of_pixels += 2;
				}
			}
		}
		// Last plane
		for (line_counter=0; line_counter < logical_y_dimension; line_counter++)
		{
			for (pixel_counter=0; pixel_counter < logical_x_dimension; pixel_counter++)
			{
				sum += real_values[address];
				address++;
			}
			address += padding_jump_value;
		}
		number_of_pixels += logical_x_dimension * logical_y_dimension;
	}
	return sum/float(number_of_pixels);
}

float Image::ReturnAverageOfRealValuesAtRadius(float wanted_mask_radius)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	double sum = 0.0;
	long address = 0;
	long number_of_pixels = 0;
	int		i;
	int		j;
	int 	k;
	float	x;
	float	y;
	float	z;
	float   mask_radius_squared;
	float	distance_from_center_squared;

	mask_radius_squared = powf(wanted_mask_radius, 2);
	number_of_pixels = 0;
	for (k = 0; k < logical_z_dimension; k++)
	{
		z = powf(k - physical_address_of_box_center_z, 2);

		for (j = 0; j < logical_y_dimension; j++)
		{
			y = powf(j - physical_address_of_box_center_y, 2);

			for (i = 0; i < logical_x_dimension; i++)
			{
				x = powf(i - physical_address_of_box_center_x, 2);

				distance_from_center_squared = x + y + z;

				if (fabsf(distance_from_center_squared -mask_radius_squared) < 4.0)
				{
					sum += real_values[address];
					number_of_pixels++;
				}
				address++;
			}
			address += padding_jump_value;
		}
	}
	if (number_of_pixels > 0)
	{
		return float(sum / number_of_pixels);
	}
	else
	{
		return 0.0;
	}
}

// Find the largest voxel value, only considering voxels which are at least a certain distance from the center and from the edge in each dimension
float Image::ReturnMaximumValue(float minimum_distance_from_center, float minimum_distance_from_edge)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i,j,k;
	int i_dist_from_center, j_dist_from_center, k_dist_from_center;
	float maximum_value = - std::numeric_limits<float>::max();
	const int last_acceptable_address_x = logical_x_dimension - minimum_distance_from_edge - 1;
	const int last_acceptable_address_y = logical_y_dimension - minimum_distance_from_edge - 1;
	const int last_acceptable_address_z = logical_z_dimension - minimum_distance_from_edge - 1;
	long address = 0;


	for (k=0;k<logical_z_dimension;k++)
	{
		if (logical_z_dimension > 1)
		{
			k_dist_from_center = abs(k - physical_address_of_box_center_z);
			if (k_dist_from_center < minimum_distance_from_center || k < minimum_distance_from_edge || k > last_acceptable_address_z)
			{
				address += logical_y_dimension * (logical_x_dimension + padding_jump_value);
				continue;
			}
		}
		for (j=0;j<logical_y_dimension;j++)
		{
			j_dist_from_center = abs(j - physical_address_of_box_center_y);
			if (j_dist_from_center < minimum_distance_from_center || j < minimum_distance_from_edge || j > last_acceptable_address_y)
			{
				address += logical_x_dimension + padding_jump_value;
				continue;
			}
			for (i=0;i<logical_x_dimension;i++)
			{
				i_dist_from_center = abs(i - physical_address_of_box_center_x);
				if (i_dist_from_center < minimum_distance_from_center || i < minimum_distance_from_edge || i > last_acceptable_address_x)
				{
					address++;
					continue;
				}

				maximum_value = std::max(maximum_value,real_values[address]);
				address++;
			}
			address += padding_jump_value;
		}
	}

	return maximum_value;
}

// Find the largest voxel value, only considering voxels which are at least a certain distance from the center and from the edge in each dimension
float Image::ReturnMinimumValue(float minimum_distance_from_center, float minimum_distance_from_edge)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i,j,k;
	int i_dist_from_center, j_dist_from_center, k_dist_from_center;
	float minimum_value = std::numeric_limits<float>::max();
	const int last_acceptable_address_x = logical_x_dimension - minimum_distance_from_edge - 1;
	const int last_acceptable_address_y = logical_y_dimension - minimum_distance_from_edge - 1;
	const int last_acceptable_address_z = logical_z_dimension - minimum_distance_from_edge - 1;
	long address = 0;


	for (k=0;k<logical_z_dimension;k++)
	{
		if (logical_z_dimension > 1)
		{
			k_dist_from_center = abs(k - physical_address_of_box_center_z);
			if (k_dist_from_center < minimum_distance_from_center || k < minimum_distance_from_edge || k > last_acceptable_address_z)
			{
				address += logical_y_dimension * (logical_x_dimension + padding_jump_value);
				continue;
			}
		}
		for (j=0;j<logical_y_dimension;j++)
		{
			j_dist_from_center = abs(j - physical_address_of_box_center_y);
			if (j_dist_from_center < minimum_distance_from_center || j < minimum_distance_from_edge || j > last_acceptable_address_y)
			{
				address += logical_x_dimension + padding_jump_value;
				continue;
			}
			for (i=0;i<logical_x_dimension;i++)
			{
				i_dist_from_center = abs(i - physical_address_of_box_center_x);
				if (i_dist_from_center < minimum_distance_from_center || i < minimum_distance_from_edge || i > last_acceptable_address_x)
				{
					address++;
					continue;
				}

				minimum_value = std::min(minimum_value,real_values[address]);
				address++;
			}
			address += padding_jump_value;
		}
	}

	return minimum_value;
}

float Image::ReturnMedianOfRealValues()
{
	long number_of_voxels = logical_x_dimension * logical_y_dimension * logical_z_dimension;
	float *buffer_array = new float[number_of_voxels];

	float median_value;

	long address = 0;
	long buffer_counter = 0;

	int		i;
	int		j;
	int 	k;

	for (k = 0; k < logical_z_dimension; k++)
	{
		for (j = 0; j < logical_y_dimension; j++)
		{
			for (i = 0; i < logical_x_dimension; i++)
			{
				buffer_array[buffer_counter] = real_values[address];

				buffer_counter++;
				address++;
			}

			address += padding_jump_value;
		}
	}

	std::sort(buffer_array, buffer_array + number_of_voxels -1);
	median_value = buffer_array[number_of_voxels / 2];
	delete [] buffer_array;

	return median_value;
}

//TODO: consolidate (reduce code duplication) by using an Empirical distribution object
float Image::ReturnAverageOfRealValues(float wanted_mask_radius, bool invert_mask)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	double sum = 0.0;
	long address = 0;
	long number_of_pixels = 0;
	int		i;
	int		j;
	int 	k;
	float	x;
	float	y;
	float	z;
	float   mask_radius_squared;
	float	distance_from_center_squared;

	if (wanted_mask_radius > 0.0)
	{
		mask_radius_squared = powf(wanted_mask_radius, 2);
		number_of_pixels = 0;
		for (k = 0; k < logical_z_dimension; k++)
		{
			z = powf(k - physical_address_of_box_center_z, 2);

			for (j = 0; j < logical_y_dimension; j++)
			{
				y = powf(j - physical_address_of_box_center_y, 2);

				for (i = 0; i < logical_x_dimension; i++)
				{
					x = powf(i - physical_address_of_box_center_x, 2);

					distance_from_center_squared = x + y + z;

					if (invert_mask)
					{
						if (distance_from_center_squared > mask_radius_squared)
						{
							sum += real_values[address];
							number_of_pixels++;
						}
					}
					else
					{
						if (distance_from_center_squared <= mask_radius_squared)
						{
							sum += real_values[address];
							number_of_pixels++;
						}
					}
					address++;
				}
				address += padding_jump_value;
			}
		}
		if (number_of_pixels > 0)
		{
			return float(sum / number_of_pixels);
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		for (k = 0; k < logical_z_dimension; k++)
		{
			for (j = 0; j < logical_y_dimension; j++)
			{
				for (i = 0; i < logical_x_dimension; i++)
				{
					sum += real_values[address];
					address++;
				}
				address += padding_jump_value;
			}
		}

	}
	return float(sum / (long(logical_x_dimension) * long(logical_y_dimension) * long(logical_z_dimension)));
}


void Image::ComputeAverageAndSigmaOfValuesInSpectrum(float minimum_radius, float maximum_radius, float &average, float &sigma, int cross_half_width)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");
	MyDebugAssertTrue(maximum_radius > minimum_radius,"Maximum radius must be greater than minimum radius");
	MyDebugAssertTrue(logical_z_dimension == 1, "Meant for images, not volumes");

	// Private variables
	int i, j;
	float x_sq, y_sq, rad_sq;
	EmpiricalDistribution my_distribution;
	const float min_rad_sq = powf(minimum_radius,2);
	const float max_rad_sq = powf(maximum_radius,2);
	const float cross_half_width_sq = powf(cross_half_width,2);
	long address = -1;

	for (j=0;j<logical_y_dimension;j++)
	{
		y_sq = powf(j-physical_address_of_box_center_y,2);
		if (y_sq <= cross_half_width_sq)
		{
			address += logical_x_dimension + padding_jump_value;
			continue;
		}
		for (i=0;i<logical_x_dimension;i++)
		{
			address++;
			x_sq = powf(i-physical_address_of_box_center_x,2);
			if (x_sq <= cross_half_width_sq) continue;
			rad_sq = x_sq + y_sq;
			if (rad_sq > min_rad_sq && rad_sq < max_rad_sq)
			{
				my_distribution.AddSampleValue(real_values[address]);
			}

		}
		address += padding_jump_value;
	}
	average = my_distribution.GetSampleMean();
	sigma = sqrtf(my_distribution.GetSampleVariance());

}

void Image::ZeroCentralPixel()
{
	MyDebugAssertTrue(logical_z_dimension == 1, "Meant for images, not volumes");

	if (is_in_real_space == false)
	{
		complex_values[0] = 0.0f * I + 0.0f;
	}
	else
	{
		int i,j;
		long address = 0;

		for (j=0;j<logical_y_dimension;j++)
		{
			for (i=0;i<logical_x_dimension;i++)
			{
				if (j==physical_address_of_box_center_y && i==physical_address_of_box_center_x)
				{
					real_values[address] = 0;
				}
				address++;
			}
			address += padding_jump_value;
		}

	}
}


void Image::SetMaximumValueOnCentralCross(float maximum_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");
	MyDebugAssertTrue(logical_z_dimension == 1, "Meant for images, not volumes");

	int i,j;
	long address = 0;

	for (j=0;j<logical_y_dimension;j++)
	{
		for (i=0;i<logical_x_dimension;i++)
		{
			if (j==physical_address_of_box_center_y || i==physical_address_of_box_center_x)
			{
				real_values[address] = std::min(maximum_value,real_values[address]);
			}
			address++;
		}
		address += padding_jump_value;
	}

}

// The image is assumed to be an amplitude spectrum, which we want to correlate with a set of CTF parameters
float Image::GetCorrelationWithCTF(CTF ctf)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");
	MyDebugAssertTrue(logical_z_dimension == 1, "Meant for images, not volumes");
	MyDebugAssertTrue(ctf.GetLowestFrequencyForFitting() > 0, "Will not work with lowest frequency for fitting of 0.");

	// Local variables
	double			cross_product = 0.0;
	double			norm_image = 0.0;
	double			norm_ctf = 0.0;
	long			number_of_values = 0;
	int				i,j;
	float				i_logi, j_logi;
	float			i_logi_sq, j_logi_sq;
	const float		inverse_logical_x_dimension = 1.0 / float(logical_x_dimension);
	const float		inverse_logical_y_dimension = 1.0 / float(logical_y_dimension);
	float			current_spatial_frequency_squared;
	const float		lowest_freq = powf(ctf.GetLowestFrequencyForFitting(),2);
	const float		highest_freq = powf(ctf.GetHighestFrequencyForFitting(),2);
	long			address = 0;
	float			current_azimuth;
	float			current_ctf_value;
	const int		central_cross_half_width = 10;
	float			astigmatism_penalty;
		
	// Loop over half of the image (ignore Friedel mates)
	for (j=0;j<logical_y_dimension;j++)
	{
		// DNM: Moved test on j out of inner loop, loop i only as far as needed, use an address computed for each line (speeds up ~13%)
		if (j < physical_address_of_box_center_y - central_cross_half_width || j > physical_address_of_box_center_y + central_cross_half_width)
		{
			address = j * (padding_jump_value + 2 * physical_address_of_box_center_x);
			j_logi = float(j-physical_address_of_box_center_y)*inverse_logical_y_dimension;
			j_logi_sq = powf(j_logi,2);
			for (i=0;i < physical_address_of_box_center_x - central_cross_half_width;i++)
			{
				i_logi = float(i-physical_address_of_box_center_x)*inverse_logical_x_dimension;
				i_logi_sq = powf(i_logi,2);
				
				// Where are we?
				current_spatial_frequency_squared = j_logi_sq + i_logi_sq;
				
				if (current_spatial_frequency_squared > lowest_freq && current_spatial_frequency_squared < highest_freq)
				{
					current_azimuth = atan2f(j_logi,i_logi);
					current_ctf_value = fabsf(ctf.Evaluate(current_spatial_frequency_squared,current_azimuth));
					// accumulate results
					number_of_values++;
					cross_product += real_values[address + i] * current_ctf_value;
					norm_image    += pow(real_values[address + i],2);
					norm_ctf      += pow(current_ctf_value,2);

				} // end of test whether within min,max frequency range
									
				// We're going to the next pixel
				//address++;
			}
			// We're going to the next line
			//address += padding_jump_value + physical_address_of_box_center_x;
		}
	}
	
	// Compute the penalty due to astigmatism
	if (ctf.GetAstigmatismTolerance() > 0.0)
	{
		astigmatism_penalty = powf(ctf.GetAstigmatism(),2) * 0.5 / powf(ctf.GetAstigmatismTolerance(),2) / float(number_of_values);
	}
	else
	{
		astigmatism_penalty = 0.0;
	}

	// The final score
	return cross_product / sqrt(norm_image * norm_ctf) - astigmatism_penalty;
}

// DNM: Set up for doing correlations with a CTF by making tables of all the pixels that are included, and their frequencies and azimuths
// When called with addresses NULL, it simply returns the number of values needed in the arrays
// When called with proper addresses, it fills the array and computes norm_image and image_mean, which are constant
void Image::SetupQuickCorrelationWithCTF(CTF ctf, int &number_of_values, double &norm_image, double &image_mean, int *addresses, float *spatial_frequency_squared, float *azimuth)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");
	MyDebugAssertTrue(logical_z_dimension == 1, "Meant for images, not volumes");
	MyDebugAssertTrue(ctf.GetLowestFrequencyForFitting() > 0, "Will not work with lowest frequency for fitting of 0.");

	// Local variables
	int				i,j;
	float				i_logi, j_logi;
	float			i_logi_sq, j_logi_sq;
	const float		inverse_logical_x_dimension = 1.0 / float(logical_x_dimension);
	const float		inverse_logical_y_dimension = 1.0 / float(logical_y_dimension);
	float			current_spatial_frequency_squared;
	
	const float		lowest_freq = powf(ctf.GetLowestFrequencyForFitting(),2);
	const float		highest_freq = powf(ctf.GetHighestFrequencyForFitting(),2);
	int			address = 0;
	float			current_azimuth;
	const int		central_cross_half_width = 10;
	double image_sum = 0.;

	number_of_values = 0;
	norm_image = 0;
	image_mean = 0.;
		
	// Loop over half of the image (ignore Friedel mates)
	for (j=0;j<logical_y_dimension;j++)
	{
		if (j < physical_address_of_box_center_y - central_cross_half_width || j > physical_address_of_box_center_y + central_cross_half_width)
		{
			address = j * (padding_jump_value + 2 * physical_address_of_box_center_x);
			j_logi = float(j-physical_address_of_box_center_y)*inverse_logical_y_dimension;
			j_logi_sq = powf(j_logi,2);
			for (i=0;i < physical_address_of_box_center_x - central_cross_half_width;i++)
			{
				i_logi = float(i-physical_address_of_box_center_x)*inverse_logical_x_dimension;
				i_logi_sq = powf(i_logi,2);
					
				// Where are we?
				current_spatial_frequency_squared = j_logi_sq + i_logi_sq;
					
				if (current_spatial_frequency_squared > lowest_freq && current_spatial_frequency_squared < highest_freq)
				{
					current_azimuth = atan2f(j_logi,i_logi);
					if (addresses)
					{
						addresses[number_of_values] = address + i;
						spatial_frequency_squared[number_of_values] = current_spatial_frequency_squared;
						azimuth[number_of_values] = current_azimuth;
						image_sum += real_values[address + i];
					}
					number_of_values++;
				} // end of test whether within min,max frequency range
					
			}
		}
	}

	// Now get sum of squared deviations from mean, more accurate than using raw cross-products
	if (addresses) 
  {
		image_mean = image_sum / number_of_values;
		for (i = 0; i < number_of_values; i++)
			norm_image += pow(real_values[addresses[i]] - image_mean, 2);
	}
}

// DNM: Computes correlation with the current CTF estimate given the pixel indexes, frequency and azimuth values
// It is about 30% faster than original and now returns true correlation coefficient
float Image::QuickCorrelationWithCTF(CTF ctf, int number_of_values, double norm_image, double image_mean, int *addresses, float *spatial_frequency_squared, float *azimuth)
{

	// Local variables
	int				i,j;
	double			cross_product = 0.0;
	double			norm_ctf = 0.0;
	double		ctf_sum = 0.;
	float			current_ctf_value;
	float			astigmatism_penalty;

	for (i = 0; i < number_of_values; i++) {
		j = addresses[i];
		current_ctf_value = fabsf(-sin(ctf.PhaseShiftGivenSquaredSpatialFrequencyAndAzimuth(spatial_frequency_squared[i], azimuth[i])));
		cross_product += real_values[j] * current_ctf_value;
		norm_ctf			+= pow(current_ctf_value,2);
		ctf_sum				+= current_ctf_value;
	}

	// Compute the penalty due to astigmatism
	if (ctf.GetAstigmatismTolerance() > 0.0)
	{
		astigmatism_penalty = powf(ctf.GetAstigmatism(),2) * 0.5 / powf(ctf.GetAstigmatismTolerance(),2) / float(number_of_values);
	}
	else
	{
		astigmatism_penalty = 0.0;
	}

	// The final score: norm_image is already a sum of squared deviations from mean; norm_ctf requires adjustment to give true CC
	return (cross_product - image_mean * ctf_sum) / sqrt(norm_image * (norm_ctf - ctf_sum * ctf_sum / number_of_values)) - astigmatism_penalty;
}

void Image::ApplyMirrorAlongY()
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");
	MyDebugAssertTrue(logical_z_dimension == 1, "Meant for images, not volumes");

	int i,j;
	long address = logical_x_dimension + padding_jump_value;
	int j_dist;
	float temp_value;

	for (j = 1; j < physical_address_of_box_center_y; j++)
	{
		j_dist = 2 * ( physical_address_of_box_center_y - j ) * (logical_x_dimension + padding_jump_value);

		for (i = 0; i < logical_x_dimension; i++)
		{
			temp_value = real_values[address];
			real_values[address] = real_values[address+j_dist];
			real_values[address+j_dist] = temp_value;
			address++;
		}
		address += padding_jump_value;
	}

	// The column j=0 is undefined, we set it to the average of the values that were there before the mirror operation was applied
	temp_value = 0;
	for (i=0 ; i < logical_x_dimension; i++) {
		temp_value += real_values[i];
	}
	temp_value /= float(logical_x_dimension);
	for (i=0 ; i < logical_x_dimension; i++) {
		real_values[i] = temp_value;
	}
}

void Image::AddImage(Image *other_image)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] += other_image->real_values[pixel_counter];
	}

}

void Image::SubtractImage(Image *other_image)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(HasSameDimensionsAs(other_image),"Images should have same dimensions, but they don't: %i %i %i        %i %i %i",logical_x_dimension,logical_y_dimension,logical_z_dimension,other_image->logical_x_dimension,other_image->logical_y_dimension,other_image->logical_z_dimension);

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] -= other_image->real_values[pixel_counter];
	}

}

void Image::SubtractSquaredImage(Image *other_image)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");

	for (long pixel_counter = 0; pixel_counter < real_memory_allocated; pixel_counter++)
	{
		real_values[pixel_counter] -= powf(other_image->real_values[pixel_counter],2);
	}

}

int Image::ReturnFourierLogicalCoordGivenPhysicalCoord_X(int physical_index)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(physical_index <= physical_upper_bound_complex_x, "index out of bounds");

    //if (physical_index >= physical_index_of_first_negative_frequency_x)
    if (physical_index > physical_address_of_box_center_x)
    {
    	 return physical_index - logical_x_dimension;
    }
    else return physical_index;
}


int Image::ReturnFourierLogicalCoordGivenPhysicalCoord_Y(int physical_index)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(physical_index <= physical_upper_bound_complex_y, "index out of bounds");

    if (physical_index >= physical_index_of_first_negative_frequency_y)
    {
    	 return physical_index - logical_y_dimension;
    }
    else return physical_index;
}

int Image::ReturnFourierLogicalCoordGivenPhysicalCoord_Z(int physical_index)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(physical_index <= physical_upper_bound_complex_z, "index out of bounds");

    if (physical_index >= physical_index_of_first_negative_frequency_z)
    {
    	 return physical_index - logical_z_dimension;
    }
    else return physical_index;
}



//  \brief  Compute the 1D rotational average
//			The first element will be the value at the center/origin of the image.
//			It is assumed the X axis of the Curve object has been setup already. It should run from 0.0 to the maximum value
//			possible, which is approximately sqrt(2)*0.5 in Fourier space or sqrt(2)*0.5*logical_dimension in real space
//			(to compute this properly, use ReturnMaximumDiagonalRadius * fourier_voxel_size). To use
//			The Fourier space radius convention in real space, give fractional_radius_in_real_space
void Image::Compute1DRotationalAverage(Curve &average, Curve &number_of_values, bool fractional_radius_in_real_space)
{

	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(average.number_of_points == number_of_values.number_of_points,"Curves do not have the same number of points");

	int i;
	int j;
	int k;
	float rad;
	long address;

	// Initialise
	average.ZeroYData();
	number_of_values.ZeroYData();
	address = 0;


	//
	if (is_in_real_space && !fractional_radius_in_real_space)
	{
		int i_logi,j_logi,k_logi;

		for (k=0;k<logical_z_dimension;k++)
		{ 
      /* IMOD icl 11 switch to multiply 3 times 
      k_logi = pow((k-physical_address_of_box_center_z),2); */
      k_logi = (k-physical_address_of_box_center_z) * (k-physical_address_of_box_center_z);
			for (j=0;j<logical_y_dimension;j++)
			{
				/*j_logi = pow((j-physical_address_of_box_center_y),2) + k_logi;*/
				j_logi = (j-physical_address_of_box_center_y) * (j-physical_address_of_box_center_y) + k_logi;
				for (i=0;i<logical_x_dimension;i++)
				{
					/*i_logi = pow((i-physical_address_of_box_center_x),2) + j_logi;*/
					i_logi = (i-physical_address_of_box_center_x) * (i-physical_address_of_box_center_x) + j_logi;
					//
					rad = sqrt(float(i_logi));
					//
					average.AddValueAtXUsingLinearInterpolation(rad,real_values[address],true);
					number_of_values.AddValueAtXUsingLinearInterpolation(rad,1.0,true);

					// Increment the address
					address ++;
				}
				// End of the line in real space
				address += padding_jump_value;
			}
		}
	}
	else
	{
		float i_logi,j_logi,k_logi;

		if (is_in_real_space && fractional_radius_in_real_space)
		{
			for (k=0;k<logical_z_dimension;k++)
			{
				k_logi = pow((k-physical_address_of_box_center_z) * fourier_voxel_size_z,2);
				for (j=0;j<logical_y_dimension;j++)
				{
					j_logi = pow((j-physical_address_of_box_center_y) * fourier_voxel_size_y,2) + k_logi;
					for (i=0;i<logical_x_dimension;i++)
					{
						i_logi = pow((i-physical_address_of_box_center_x) * fourier_voxel_size_x,2) + j_logi;
						//
						rad = sqrt(float(i_logi));
						//
						average.AddValueAtXUsingLinearInterpolation(rad,real_values[address],true);
						number_of_values.AddValueAtXUsingLinearInterpolation(rad,1.0,true);

						// Increment the address
						address ++;
					}
					// End of the line in real space
					address += padding_jump_value;
				}
			}
		}
		else
		{
			for (k=0;k<logical_z_dimension;k++)
			{
				k_logi = pow(ReturnFourierLogicalCoordGivenPhysicalCoord_Z(k) * fourier_voxel_size_z,2);
				for (j=0;j<logical_y_dimension;j++)
				{
					j_logi = pow(ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * fourier_voxel_size_y,2) + k_logi;
					for (i=0;i<physical_upper_bound_complex_x;i++)
					{
						i_logi = pow(i * fourier_voxel_size_x,2) + j_logi;
						//
						if (FourierComponentIsExplicitHermitianMate(i,j,k)) continue;
						rad = sqrt(float(i_logi));
						//
						average.AddValueAtXUsingLinearInterpolation(rad,abs(complex_values[address]),true);
						number_of_values.AddValueAtXUsingLinearInterpolation(rad,1.0,true);

						// Increment the address
						address ++;
					}
				}
			}
		}
	}

	// Do the actual averaging
	for (int counter = 0; counter < average.number_of_points; counter ++ )
	{
		if (number_of_values.data_y[counter] != 0.0) average.data_y[counter] /=number_of_values.data_y[counter];
	}
}

// It is assumed the curve objects are already setup with an X axis in reciprocal pixels (i.e. origin is 0.0, Nyquist is 0.5)
void Image::Compute1DPowerSpectrumCurve(Curve *curve_with_average_power, Curve *curve_with_number_of_values)
{

	MyDebugAssertTrue(is_in_memory,"Memory not allocated");
	MyDebugAssertFalse(is_in_real_space,"Image not in Fourier space");
	MyDebugAssertTrue(curve_with_average_power->number_of_points > 0, "Curve not setup");
	MyDebugAssertTrue(curve_with_average_power->data_x[0] == 0.0, "Curve does not start at x = 0\n");
	MyDebugAssertTrue(curve_with_average_power->data_x[curve_with_average_power->number_of_points-1] >= 0.5, "Curve does not go to at least x = 0.5 (it goes to %f)\n",curve_with_average_power->data_x[curve_with_average_power->number_of_points-1]);
	MyDebugAssertTrue(curve_with_average_power->number_of_points == curve_with_number_of_values->number_of_points, "Curves need to have the same number of points");
	MyDebugAssertTrue(curve_with_average_power->data_x[0] == curve_with_number_of_values->data_x[0], "Curves need to have the same starting point");
	MyDebugAssertTrue(curve_with_average_power->data_x[curve_with_average_power->number_of_points-1] == curve_with_number_of_values->data_x[curve_with_number_of_values->number_of_points-1], "Curves need to have the same ending point");


	int i,j,k;
	float sq_dist_x, sq_dist_y, sq_dist_z;
	int counter;
	long address;
	float spatial_frequency;
	int number_of_hermitian_mates = 0;

	// Make sure the curves are clean
	curve_with_average_power->ZeroYData();
	curve_with_number_of_values->ZeroYData();


	// Get amplitudes and sum them into the curve object
	address = 0;
	for ( k = 0; k <= physical_upper_bound_complex_z; k ++ )
	{
		sq_dist_z = powf(ReturnFourierLogicalCoordGivenPhysicalCoord_Z(k) * fourier_voxel_size_z,2);
		for ( j = 0; j <= physical_upper_bound_complex_y; j ++ )
		{
			sq_dist_y = powf(ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j) * fourier_voxel_size_y,2);
			for  ( i = 0; i <= physical_upper_bound_complex_x; i ++ )
			{
				if (FourierComponentIsExplicitHermitianMate(i,j,k))
				{
					number_of_hermitian_mates++;
					address ++;
					continue;
				}
				else
				{
					sq_dist_x = powf(i * fourier_voxel_size_x,2);
					spatial_frequency = sqrtf(sq_dist_x+sq_dist_y+sq_dist_z);

					// TODO: this could be made faster by doing both interpolations in one go, so one wouldn't have to work out twice between which points the interpolation will happen
					curve_with_average_power->AddValueAtXUsingLinearInterpolation(spatial_frequency,real(complex_values[address]) * real(complex_values[address]) + imag(complex_values[address]) * imag(complex_values[address]), true );
					curve_with_number_of_values->AddValueAtXUsingLinearInterpolation(spatial_frequency,1.0, true);

					address ++;
				}
			}
		}
	}

	// Do the actual averaging
	for ( counter = 0; counter < curve_with_average_power->number_of_points; counter ++ )
	{
		if ( curve_with_number_of_values->data_y[counter] > 0.0 )
		{
			curve_with_average_power->data_y[counter] /= curve_with_number_of_values->data_y[counter];
		}
		else
		{
			curve_with_average_power->data_y[counter] = 0.0;
		}
	}

}


void Image::ComputeAmplitudeSpectrumFull2D(Image *amplitude_spectrum)
{
	MyDebugAssertTrue(is_in_memory,"Memory not allocated");
	MyDebugAssertTrue(amplitude_spectrum->is_in_memory, "Other image memory not allocated");
	MyDebugAssertTrue(HasSameDimensionsAs(amplitude_spectrum), "Images do not have same dimensions");
	MyDebugAssertFalse(is_in_real_space,"Image not in Fourier space");

	int ampl_addr_i;
	int ampl_addr_j;
	int image_addr_i;
	int image_addr_j;
	int i_mate;
	int j_mate;

	long address_in_amplitude_spectrum = 0;
	long address_in_self;

	// Loop over the amplitude spectrum
	for (ampl_addr_j = 0; ampl_addr_j < amplitude_spectrum->logical_y_dimension; ampl_addr_j++)
	{
		for (ampl_addr_i = 0; ampl_addr_i < amplitude_spectrum->logical_x_dimension; ampl_addr_i++)
		{
			address_in_self = ReturnFourier1DAddressFromLogicalCoord(ampl_addr_i-amplitude_spectrum->physical_address_of_box_center_x,ampl_addr_j-amplitude_spectrum->physical_address_of_box_center_y,0);
			amplitude_spectrum->real_values[address_in_amplitude_spectrum] = abs(complex_values[address_in_self]);
			address_in_amplitude_spectrum++;
		}
		address_in_amplitude_spectrum += amplitude_spectrum->padding_jump_value;
	}

	// Done
	amplitude_spectrum->is_in_real_space = true;
	amplitude_spectrum->object_is_centred_in_box = true;
}


/*
 * Real-space box convolution meant for 2D amplitude spectra
 *
 * This is adapted from the MSMOOTH subroutine from CTFFIND3, with a different wrap-around behaviour.
 * Also, in this version, we loop over the full 2D, rather than just half - this runs faster because less logic within the loop
 * DNM rewrote this to be vastly faster
 */
void Image::SpectrumBoxConvolution(Image *output_image, int box_size, float minimum_radius)
{
	MyDebugAssertTrue(IsEven(box_size) == false,"Box size must be odd");
	MyDebugAssertTrue(logical_z_dimension == 1,"Volumes not supported");
	MyDebugAssertTrue(output_image->is_in_memory == true,"Output image not allocated");
	MyDebugAssertTrue(HasSameDimensionsAs(output_image),"Output image does not have same dimensions as image");

	// Variables
	const int half_box_size = (box_size-1)/2;
	const int cross_half_width_to_ignore = 1;
	int i;
	int i_sq;
	int ii;
	int j;
	int j_sq;
	int jj;
	int num_voxels;
	int m;
	int l;
	const float minimum_radius_sq = pow(minimum_radius,2);
	float radius_sq;
	const int first_i_to_ignore = physical_address_of_box_center_x - cross_half_width_to_ignore;
	const int last_i_to_ignore  = physical_address_of_box_center_x + cross_half_width_to_ignore;
	const int first_j_to_ignore = physical_address_of_box_center_y - cross_half_width_to_ignore;
	const int last_j_to_ignore  = physical_address_of_box_center_y + cross_half_width_to_ignore;

	// Addresses
	long address_within_output = 0;
	long address_within_input;
 
  // Starting and ending x indexes of one or two loops for each line
  int *x1start = new int[logical_x_dimension];
  int *x1end = new int[logical_x_dimension];
  int *x2start = new int[logical_x_dimension];
  int *x2end = new int[logical_x_dimension];
  int *numInLineSum = new int[logical_x_dimension];
  float *lineSums = new float[logical_x_dimension * logical_y_dimension];
  float sum;
  int ybase;

  // Get the limits for one or two loops for making line sums at each X position
  for (i = 0; i < logical_x_dimension; i++) 
  {
    x1start[i] = i - half_box_size;
    x1end[i] = i + half_box_size;
    x2start[i] = 0;
    x2end[i] = -1;

    // Wrap around left edge
    if (x1start[i] < 0) 
    {
      x2start[i] = x1start[i] + logical_x_dimension;
      x2end[i] = logical_x_dimension - 1;
      x1start[i] = 0;
    } 

    // Or wrap around right edge
    else if (x1end[i] >= logical_x_dimension)
    {
      x2end[i] = x1end[i] - logical_x_dimension;
      x2start[i] = 0;
      x1end[i] = logical_x_dimension - 1;
    }

    // Or handle intersection with the central cross by trimming or splitting into two loops
    else if (x1start[i] <= last_i_to_ignore && x1end[i] >= first_i_to_ignore)
    {
      if (x1start[i] >= first_i_to_ignore) 
        x1start[i] = last_i_to_ignore + 1;
      else if (x1end[i] <= last_i_to_ignore)
        x1end[i] = first_i_to_ignore - 1;
      else
      {
        x2end[i] = x1end[i];
        x2start[i] = last_i_to_ignore + 1;
        x1end[i] = first_i_to_ignore - 1;
      }
    }
    numInLineSum[i] = x1end[i] + 1 - x1start[i];
    if (x2end[i] >= x2start[i])
      numInLineSum[i] += x2end[i] + 1 - x2start[i];
  }

  // Loop over Y positions for line sums
	for (jj = 0; jj < logical_y_dimension; jj++)
  {
    ybase = jj * (logical_x_dimension + padding_jump_value);

    // Form line sums at each X position
		for (i = 0; i < logical_x_dimension; i++)
    {
      sum = 0.;
      for (ii = x1start[i]; ii <= x1end[i]; ii++)
        sum += real_values[ii + ybase];
      for (ii = x2start[i]; ii <= x2end[i]; ii++)
        sum += real_values[ii + ybase];
      lineSums[i + jj * logical_x_dimension] = sum;
    }
  }

	// Loop over the output image
	for (j = 0; j < logical_y_dimension; j++)
	{
    /* IMOD icl 11 switch pow to multiply twice
    j_sq = pow((j - physical_address_of_box_center_y),2);*/
		j_sq = (j - physical_address_of_box_center_y) * (j - physical_address_of_box_center_y);

		for (i = 0; i < logical_x_dimension; i++)
		{
      /*i_sq = pow((i - physical_address_of_box_center_x),2); */
			i_sq = (i - physical_address_of_box_center_x) * (i - physical_address_of_box_center_x);

			radius_sq = float(i_sq+j_sq);

			if ( radius_sq <= minimum_radius_sq )
			{
				output_image->real_values[address_within_output] = real_values[address_within_output];
			}
			else
			{
				output_image->real_values[address_within_output] = 0.0e0;
				num_voxels = 0;

        // Loop over the lines to sum at this pixel to get the box sum
				for ( m = - half_box_size; m <= half_box_size; m++)
				{
					jj = j + m;
					// wrap around
					if (jj < 0) { jj += logical_y_dimension; }
					if (jj >= logical_y_dimension) { jj -= logical_y_dimension; }

					// In central cross?
					//if ( abs(jj - physical_address_of_box_center_y) <= cross_half_width_to_ignore ) { continue; }
					if ( jj >= first_j_to_ignore && jj <= last_j_to_ignore) { continue; }

          output_image->real_values[address_within_output] += lineSums[i + jj * logical_x_dimension];
          num_voxels += numInLineSum[i];

				} // end of loop over the box

				if (num_voxels == 0)
				{
          // DNM: if it happens, surely that should be from same address not whatever address_within_input was
					output_image->real_values[address_within_output] = real_values[address_within_output];
				}
				else
				{
					output_image->real_values[address_within_output] /= float(num_voxels);
				}
			}

			address_within_output++;
		}
		address_within_output += output_image->padding_jump_value;
	}

  delete [] x1start;
  delete [] x1end;
  delete [] x2start;
  delete [] x2end;
  delete [] lineSums;
  delete [] numInLineSum;
}



/*

void Image::SpectrumBoxConvolution(Image *output_image, int box_size, float minimum_radius)
{
	MyDebugAssertTrue(IsEven(box_size) == false,"Box size must be odd");
	MyDebugAssertTrue(logical_z_dimension == 1,"Volumes not supported");
	MyDebugAssertTrue(output_image->is_in_memory == true,"Output image not allocated");
	MyDebugAssertTrue(HasSameDimensionsAs(output_image),"Output image does not have same dimensions as image");

	// Variables
	int half_box_size = (box_size-1)/2;
	int cross_half_width_to_ignore = 1;
	int i;
	int i_friedel;
	int i_sq;
	int ii;
	int ii_friedel;
	int ii_sq;
	int iii;
	int j;
	int j_friedel;
	int j_sq;
	int jj;
	int jj_friedel;
	int jj_sq;
	int jjj;
	float radius;
	int num_voxels;
	int m;
	int l;

	// Addresses
	long address_within_output = 0;
	long address_within_input;

	// Loop over the output image. To save time, we only loop over one half of the image [BUG: actually this is looping over the full image!
	for (j = 0; j < logical_y_dimension; j++)
	{
		j_friedel = 2 * physical_address_of_box_center_y - j;
		j_sq = powf((j - physical_address_of_box_center_y),2);

		for (i = 0; i < logical_x_dimension; i++)
		{
			i_friedel = 2 * physical_address_of_box_center_x - i;
			i_sq = powf((i - physical_address_of_box_center_x),2);

			//address_within_output = ReturnReal1DAddressFromPhysicalCoord(i,j,0);

			radius = sqrt(float(i_sq+j_sq));

			if ( radius <= minimum_radius )
			{
				output_image->real_values[address_within_output] = real_values[address_within_output];
			}
			else
			{
				output_image->real_values[address_within_output] = 0.0e0;
				num_voxels = 0;

				for ( m = - half_box_size; m <= half_box_size; m++)
				{
					jj = j + m;
					if (jj < 0) { jj += logical_y_dimension; }
					if (jj >= logical_y_dimension) { jj -= logical_y_dimension; }
					jj_friedel = 2 * physical_address_of_box_center_y - jj;
					jj_sq = powf((jj - physical_address_of_box_center_y),2);

					for ( l = - half_box_size; l <= half_box_size; l++)
					{
						ii = i + l;
						if (ii < 0) { ii += logical_x_dimension; }
						if (ii >= logical_x_dimension) { ii -= logical_x_dimension; }
						ii_friedel = 2 * physical_address_of_box_center_x - ii;
						ii_sq = powf((ii - physical_address_of_box_center_x),2);

						// Friedel or not?
						if ( ii > physical_address_of_box_center_x)
						{
							iii = ii_friedel;
							jjj = jj_friedel;
							if (jjj > logical_y_dimension - 1 || iii > logical_x_dimension - 1) { continue; }
						}
						else
						{
							iii = ii;
							jjj = jj;
						}

						// In central cross?
						if ( abs(iii - physical_address_of_box_center_x) <= cross_half_width_to_ignore || abs(jjj - physical_address_of_box_center_y) <= cross_half_width_to_ignore ) { continue; }

						address_within_input = ReturnReal1DAddressFromPhysicalCoord(iii,jjj,0);

						if ( iii < logical_x_dimension && jjj < logical_y_dimension ) // it sometimes happens that we end up on Nyquist Friedel mates that we don't have (perhaps this can be fixed)
						{
							output_image->real_values[address_within_output] += real_values[address_within_input];
						}
						num_voxels++; // not sure why this is not within the if branch, like the addition itself - is this a bug?

					}
				} // end of loop over the box

				if (num_voxels == 0)
				{
					output_image->real_values[address_within_output] = real_values[address_within_input];
				}
				else
				{
					output_image->real_values[address_within_output] /= float(num_voxels);
				}
			}

			if (j_friedel < logical_y_dimension && i_friedel < logical_x_dimension)
			{
				output_image->real_values[ReturnReal1DAddressFromPhysicalCoord(i_friedel,j_friedel,0)] = output_image->real_values[ReturnReal1DAddressFromPhysicalCoord(i,j,0)];
			}
			address_within_output++;
		}
		address_within_output += output_image->padding_jump_value;
	}

	// There are a few pixels that are not set by the logical above
	for (i = physical_address_of_box_center_x + 1; i < logical_x_dimension; i++)
	{
		i_friedel = 2 * physical_address_of_box_center_x - i;
		output_image->real_values[ReturnReal1DAddressFromPhysicalCoord(i,0,0)]                     = output_image->real_values[ReturnReal1DAddressFromPhysicalCoord(i_friedel,0,0)];
		output_image->real_values[ReturnReal1DAddressFromPhysicalCoord(i,logical_y_dimension-1,0)] = output_image->real_values[ReturnReal1DAddressFromPhysicalCoord(i_friedel,logical_y_dimension - 1,0)];
	}
}
*/

//#pragma GCC push_options
//#pragma GCC optimize ("O0")
// Taper edges of image so that there are no sharp discontinuities in real space
// This is a re-implementation of the MRC program taperedgek.for (Richard Henderson, 1987)
/* IMOD: Unneeded and has an abort()
void Image::TaperEdges()
{
	MyDebugAssertTrue(is_in_memory,"Image not in memory");
	MyDebugAssertTrue(is_in_real_space,"Image not in real space");

	// Private variables
	const float				fraction_averaging = 30.0;
	const float				fraction_tapering  = 30.0;
	const int				averaging_strip_width_x	=	int(logical_x_dimension/fraction_averaging); //100
	const int				averaging_strip_width_y	=	int(logical_y_dimension/fraction_averaging);
	const int				averaging_strip_width_z =   int(logical_z_dimension/fraction_averaging);
	const int				tapering_strip_width_x	=	int(logical_x_dimension/fraction_tapering); //500
	const int				tapering_strip_width_y	=	int(logical_y_dimension/fraction_tapering);
	const int				tapering_strip_width_z	=	int(logical_z_dimension/fraction_tapering);
	const int				smoothing_half_width_x	=	1; // 1
	const int				smoothing_half_width_y	=	1;
	const int				smoothing_half_width_z	=	1;
	int						current_dimension;
	int						number_of_dimensions;
	int						second_dimension;
	int						third_dimension;
	int						logical_current_dimension;
	int						logical_second_dimension;
	int						logical_third_dimension;
	int						current_tapering_strip_width;
	int						i,j,k;
	int						j_shift,k_shift;
	int						jj,kk;
	int						number_of_values_in_running_average;
	long					address;
	int						smoothing_half_width_third_dimension;
	int						smoothing_half_width_second_dimension;
	// 2D arrays
	float					*average_for_current_edge_start = NULL;
	float					*average_for_current_edge_finish = NULL;
	float					*average_for_current_edge_average = NULL;
	float					*smooth_average_for_current_edge_start = NULL;
	float					*smooth_average_for_current_edge_finish = NULL;

	// Start work


	// Check dimensions of image are OK
	if (logical_x_dimension < 2 * tapering_strip_width_x || logical_y_dimension < 2 * tapering_strip_width_y)
	{
		MyPrintWithDetails("X,Y dimensions of input image are too small: %i %i\n", logical_x_dimension,logical_y_dimension);
		abort();
	}
	if (logical_z_dimension > 1 && logical_z_dimension < 2 * tapering_strip_width_z)
	{
		MyPrintWithDetails("Z dimension is too small: %i\n",logical_z_dimension);
		abort();
	}

	if ( logical_z_dimension > 1 )
	{
		number_of_dimensions = 3;
	}
	else
	{
		number_of_dimensions = 2;
	}


	for (current_dimension=1; current_dimension <= number_of_dimensions; current_dimension++)
	{
		switch(current_dimension)
		{
		case(1):
			second_dimension = 2;
			third_dimension = 3;
			logical_current_dimension = logical_x_dimension;
			logical_second_dimension = logical_y_dimension;
			logical_third_dimension = logical_z_dimension;
			current_tapering_strip_width = tapering_strip_width_x;
			smoothing_half_width_second_dimension = smoothing_half_width_y;
			smoothing_half_width_third_dimension = smoothing_half_width_z;
			break;
		case(2):
			second_dimension = 1;
			third_dimension = 3;
			logical_current_dimension = logical_y_dimension;
			logical_second_dimension = logical_x_dimension;
			logical_third_dimension = logical_z_dimension;
			current_tapering_strip_width = tapering_strip_width_y;
			smoothing_half_width_second_dimension = smoothing_half_width_x;
			smoothing_half_width_third_dimension = smoothing_half_width_z;
			break;
		case(3):
			second_dimension = 1;
			third_dimension = 2;
			logical_current_dimension = logical_z_dimension;
			logical_second_dimension = logical_x_dimension;
			logical_third_dimension = logical_y_dimension;
			current_tapering_strip_width = tapering_strip_width_z;
			smoothing_half_width_second_dimension = smoothing_half_width_x;
			smoothing_half_width_third_dimension = smoothing_half_width_y;
			break;
		}

		// Allocate memory
		if (average_for_current_edge_start != NULL) {
			delete [] average_for_current_edge_start;
			delete [] average_for_current_edge_finish;
			delete [] average_for_current_edge_average;
			delete [] smooth_average_for_current_edge_start;
			delete [] smooth_average_for_current_edge_finish;
		}
		average_for_current_edge_start 			= new float[logical_second_dimension*logical_third_dimension];
		average_for_current_edge_finish 		= new float[logical_second_dimension*logical_third_dimension];
		average_for_current_edge_average 		= new float[logical_second_dimension*logical_third_dimension];
		smooth_average_for_current_edge_start	= new float[logical_second_dimension*logical_third_dimension];
		smooth_average_for_current_edge_finish	= new float[logical_second_dimension*logical_third_dimension];

		// Initialise memory
		for(i=0;i<logical_second_dimension*logical_third_dimension;i++)
		{
			average_for_current_edge_start[i] = 0.0;
			average_for_current_edge_finish[i] = 0.0;
			average_for_current_edge_average[i] = 0.0;
			smooth_average_for_current_edge_start[i] = 0.0;
			smooth_average_for_current_edge_finish[i] = 0.0;
		}

		//
		// Deal with X=0 and X=logical_x_dimension edges
		//
		i = 1;
		for (k=1;k<=logical_third_dimension;k++)
		{
			for (j=1;j<=logical_second_dimension;j++)
			{
				switch(current_dimension)
				{
				case(1):
						for (i=1;i<=averaging_strip_width_x;i++)
						{
							average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension] += real_values[ReturnReal1DAddressFromPhysicalCoord(i-1,j-1,k-1)];
							average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] += real_values[ReturnReal1DAddressFromPhysicalCoord(logical_x_dimension-i,j-1,k-1)];
						}
						average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension]  /= float(averaging_strip_width_x);
						average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension]  /= float(averaging_strip_width_x);
						break;
				case(2):
						for (i=1;i<=averaging_strip_width_y;i++)
						{
							average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension] += real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,i-1,k-1)];
							average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] += real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,logical_y_dimension-i,k-1)];
						}
						average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension]  /= float(averaging_strip_width_y);
						average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension]  /= float(averaging_strip_width_y);
						break;
				case(3):
						for (i=1;i<=averaging_strip_width_z;i++)
						{
							average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension] += real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,k-1,i-1)];
							average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] += real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,k-1,logical_z_dimension-i)];
						}
						average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension]  /= float(averaging_strip_width_z);
						average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension]  /= float(averaging_strip_width_z);
						break;
				}
			}
		}

		for (address=0;address<logical_second_dimension*logical_third_dimension;address++)
		{
			average_for_current_edge_average[address] = 0.5 * ( average_for_current_edge_finish[address] + average_for_current_edge_start[address]);
			average_for_current_edge_start[address] -= average_for_current_edge_average[address];
			average_for_current_edge_finish[address] -= average_for_current_edge_average[address];
		}

		// Apply smoothing parallel to edge in the form of a running average
		for (k=1;k<=logical_third_dimension;k++)
		{
			for (j=1;j<=logical_second_dimension;j++)
			{
				number_of_values_in_running_average = 0;
				// Loop over neighbourhood of non-smooth arrays
				for (k_shift=-smoothing_half_width_third_dimension;k_shift<=smoothing_half_width_third_dimension;k_shift++)
				{
					kk = k+k_shift;
					if (kk < 1 || kk > logical_third_dimension) continue;
					for (j_shift=-smoothing_half_width_second_dimension;j_shift<=smoothing_half_width_second_dimension;j_shift++)
					{
						jj = j+j_shift;
						if (jj<1 || jj > logical_second_dimension) continue;
						number_of_values_in_running_average++;

						smooth_average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension] += average_for_current_edge_start [(jj-1)+(kk-1)*logical_second_dimension];
						smooth_average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] += average_for_current_edge_finish[(jj-1)+(kk-1)*logical_second_dimension];
					}
				}
				// Now we can compute the average
				smooth_average_for_current_edge_start [(j-1)+(k-1)*logical_second_dimension] /= float(number_of_values_in_running_average);
				smooth_average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] /= float(number_of_values_in_running_average);
			}
		}

		// Taper the image
		for (i=1;i<=logical_current_dimension;i++)
		{
			if (i<=current_tapering_strip_width)
			{
				switch(current_dimension)
				{
				case(1):
						for (k=1;k<=logical_third_dimension;k++)
						{
							for (j=1;j<=logical_second_dimension;j++)
							{
								real_values[ReturnReal1DAddressFromPhysicalCoord(i-1,j-1,k-1)] -= smooth_average_for_current_edge_start[(j-1)+(k-1)*logical_second_dimension] * float(current_tapering_strip_width-i+1) / float(current_tapering_strip_width);
							}
						}
						break;
				case(2):
						for (k=1;k<=logical_third_dimension;k++)
						{
							for (j=1;j<=logical_second_dimension;j++)
							{
								real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,i-1,k-1)] -= smooth_average_for_current_edge_start[(j-1)+(k-1)*logical_second_dimension] * float(current_tapering_strip_width-i+1) / float(current_tapering_strip_width);
							}
						}
						break;
				case(3):
						for (k=1;k<=logical_third_dimension;k++)
						{
							for (j=1;j<=logical_second_dimension;j++)
							{
								real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,k-1,i-1)] -= smooth_average_for_current_edge_start[(j-1)+(k-1)*logical_second_dimension] * float(current_tapering_strip_width-i+1) / float(current_tapering_strip_width);
							}
						}
						break;
				}
			}
			else if(i >= logical_current_dimension - current_tapering_strip_width+1)
			{
				switch(current_dimension)
				{
				case(1):
						for (k=1;k<=logical_third_dimension;k++)
						{
							for (j=1;j<=logical_second_dimension;j++)
							{
								real_values[ReturnReal1DAddressFromPhysicalCoord(i-1,j-1,k-1)] -= smooth_average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] * float(current_tapering_strip_width+i-logical_current_dimension) / float(current_tapering_strip_width);
							}
						}
						break;
				case(2):
						for (k=1;k<=logical_third_dimension;k++)
						{
							for (j=1;j<=logical_second_dimension;j++)
							{
								real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,i-1,k-1)] -= smooth_average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] * float(current_tapering_strip_width+i-logical_current_dimension) / float(current_tapering_strip_width);
							}
						}
						break;
				case(3):
						for (k=1;k<=logical_third_dimension;k++)
						{
							for (j=1;j<=logical_second_dimension;j++)
							{
								real_values[ReturnReal1DAddressFromPhysicalCoord(j-1,k-1,i-1)] -= smooth_average_for_current_edge_finish[(j-1)+(k-1)*logical_second_dimension] * float(current_tapering_strip_width+i-logical_current_dimension) / float(current_tapering_strip_width);
							}
						}
						break;
				}
			}
		}

	} // end of loop over dimensions

	// Cleanup
	delete [] average_for_current_edge_start;
	delete []average_for_current_edge_finish;
	delete [] average_for_current_edge_average;
	delete [] smooth_average_for_current_edge_start;
	delete [] smooth_average_for_current_edge_finish;


}*/
//#pragma GCC pop_options



// An alternative to ClipInto which only works for 2D real space clipping into larger image. Should be faster.
void Image::ClipIntoLargerRealSpace2D(Image *other_image, float wanted_padding_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(other_image->is_in_memory, "Other image Memory not allocated");
	MyDebugAssertTrue(is_in_real_space,"Image must be in real space");
	MyDebugAssertTrue(object_is_centred_in_box, "real space image, not centred in box");
	MyDebugAssertTrue(logical_z_dimension == 1,"Image must be 2D");
	MyDebugAssertTrue(logical_x_dimension <= other_image->logical_x_dimension && logical_y_dimension <= other_image->logical_y_dimension, "Image must be smaller than other image");

	other_image->is_in_real_space = is_in_real_space;
	other_image->object_is_centred_in_box = object_is_centred_in_box;


	// Looping variables
	long address_in_self = 0;
	long address_in_other = 0;

	int i;
	int j;

	// The address boudaries in the other_image for the input image data
	// If we are clipping a (2,2) image into a (4,4) image, we should be
	// copying into addresses 1 to 2 in both directions
	// If we are clipping a logical dimension of 2 into a dimension of 5,
	// we are copying into addresses 1 to 2
	const int i_lower_bound = other_image->physical_address_of_box_center_x - physical_address_of_box_center_x;
	const int j_lower_bound = other_image->physical_address_of_box_center_y - physical_address_of_box_center_y;
	const int i_upper_bound = i_lower_bound + logical_x_dimension - 1;
	const int j_upper_bound = j_lower_bound + logical_y_dimension - 1;

	// Loop over the other (larger) image
	for (j = 0; j < other_image->logical_y_dimension; j++)
	{
		// Check whether this line is outside of the original image
		if (j < j_lower_bound || j > j_upper_bound)
		{
			// Fill this line with the padding value
			for (i = 0; i < other_image->logical_x_dimension; i++)
			{
				other_image->real_values[address_in_other] = wanted_padding_value;
				address_in_other ++;
			}
		}
		else
		{
			// This line is within the central region
			for (i = 0; i < other_image->logical_x_dimension; i++)
			{
				if (i < i_lower_bound || i > i_upper_bound)
				{
					// We are near the beginning or the end of the line
					other_image->real_values[address_in_other] = wanted_padding_value;
				}
				else
				{
					other_image->real_values[address_in_other] = real_values[address_in_self];
					address_in_self++;
				}
				address_in_other ++;
			}
		}
		// We've reached the end of the line
		address_in_other += other_image->padding_jump_value;
		if (j >= j_lower_bound) address_in_self += padding_jump_value;
	}

}

// If you don't want to clip from the center, you can give wanted_coordinate_of_box_center_{x,y,z}. This will define the pixel in the image at which other_image will be centered. (0,0,0) means center of image.
void Image::ClipInto(Image *other_image, float wanted_padding_value, bool fill_with_noise, float wanted_noise_sigma, int wanted_coordinate_of_box_center_x, int wanted_coordinate_of_box_center_y, int wanted_coordinate_of_box_center_z)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(other_image->is_in_memory, "Other image Memory not allocated");
	MyDebugAssertFalse(is_in_real_space == true && fill_with_noise == true, "Fill with noise, only for fourier space");
	MyDebugAssertFalse((! is_in_real_space) && (wanted_coordinate_of_box_center_x != 0 || wanted_coordinate_of_box_center_y != 0 || wanted_coordinate_of_box_center_z != 0), "Cannot clip off-center in Fourier space");


	long pixel_counter = 0;
	int array_address = 0;

	int temp_logical_x;
	int temp_logical_y;
	int temp_logical_z;

	int kk;
	int k;
	int kk_logi;

	int jj;
	int jj_logi;
	int j;

	int ii;
	int ii_logi;
	int i;

	double junk;

	// take other following attributes

	other_image->is_in_real_space = is_in_real_space;
	other_image->object_is_centred_in_box = object_is_centred_in_box;

	if (is_in_real_space == true)
	{
		MyDebugAssertTrue(object_is_centred_in_box, "real space image, not centred in box");

		for (kk = 0; kk < other_image->logical_z_dimension; kk++)
		{
			kk_logi = kk - other_image->physical_address_of_box_center_z;
			k = physical_address_of_box_center_z + wanted_coordinate_of_box_center_z + kk_logi;

			for (jj = 0; jj < other_image->logical_y_dimension; jj++)
			{
				jj_logi = jj - other_image->physical_address_of_box_center_y;
				j = physical_address_of_box_center_y + wanted_coordinate_of_box_center_y + jj_logi;

				for (ii = 0; ii < other_image->logical_x_dimension; ii++)
				{
					ii_logi = ii - other_image->physical_address_of_box_center_x;
					i = physical_address_of_box_center_x + wanted_coordinate_of_box_center_x + ii_logi;

					if (k < 0 || k >= logical_z_dimension || j < 0 || j >= logical_y_dimension || i < 0 || i >= logical_x_dimension)
					{
						other_image->real_values[pixel_counter] = wanted_padding_value;
					}
					else
					{
						other_image->real_values[pixel_counter] = ReturnRealPixelFromPhysicalCoord(i, j, k);
					}

					pixel_counter++;
				}

				pixel_counter+=other_image->padding_jump_value;
			}
		}
	}
	else
	{
		for (kk = 0; kk <= other_image->physical_upper_bound_complex_z; kk++)
		{
			temp_logical_z = other_image->ReturnFourierLogicalCoordGivenPhysicalCoord_Z(kk);

			//if (temp_logical_z > logical_upper_bound_complex_z || temp_logical_z < logical_lower_bound_complex_z) continue;

			for (jj = 0; jj <= other_image->physical_upper_bound_complex_y; jj++)
			{
				temp_logical_y = other_image->ReturnFourierLogicalCoordGivenPhysicalCoord_Y(jj);

				//if (temp_logical_y > logical_upper_bound_complex_y || temp_logical_y < logical_lower_bound_complex_y) continue;

				for (ii = 0; ii <= other_image->physical_upper_bound_complex_x; ii++)
				{
					temp_logical_x = ii;

					//if (temp_logical_x > logical_upper_bound_complex_x || temp_logical_x < logical_lower_bound_complex_x) continue;

					if (fill_with_noise == false) other_image->complex_values[pixel_counter] = ReturnComplexPixelFromLogicalCoord(temp_logical_x, temp_logical_y, temp_logical_z, wanted_padding_value + I * 0.0f);
					else
					{

						if (temp_logical_x < logical_lower_bound_complex_x || temp_logical_x > logical_upper_bound_complex_x || temp_logical_y < logical_lower_bound_complex_y ||temp_logical_y > logical_upper_bound_complex_y || temp_logical_z < logical_lower_bound_complex_z || temp_logical_z > logical_upper_bound_complex_z)
						{
							other_image->complex_values[pixel_counter] = (global_random_number_generator.GetNormalRandom() * wanted_noise_sigma) + (I * global_random_number_generator.GetNormalRandom() * wanted_noise_sigma);
						}
						else
						{
							other_image->complex_values[pixel_counter] = complex_values[ReturnFourier1DAddressFromLogicalCoord(temp_logical_x,temp_logical_y, temp_logical_z)];

						}


					}
					pixel_counter++;

				}

			}
		}


		// When we are clipping into a larger volume in Fourier space, there is a half-plane (vol) or half-line (2D image) at Nyquist for which FFTW
		// does not explicitly tell us the values. We need to fill them in.
		if (logical_y_dimension < other_image->logical_y_dimension || logical_z_dimension < other_image->logical_z_dimension)
		{
			// For a 2D image
			if (logical_z_dimension == 1)
			{
				jj = physical_index_of_first_negative_frequency_y;
				for (ii = 0; ii <= physical_upper_bound_complex_x; ii++)
				{
					other_image->complex_values[other_image->ReturnFourier1DAddressFromPhysicalCoord(ii,jj,0)] = complex_values[ReturnFourier1DAddressFromPhysicalCoord(ii,jj,0)];
				}
			}
			// For a 3D volume
			else
			{

				// Deal with the positive Nyquist of the 2nd dimension
				for (kk_logi = logical_lower_bound_complex_z; kk_logi <= logical_upper_bound_complex_z; kk_logi ++)
				{
					jj = physical_index_of_first_negative_frequency_y;
					jj_logi = logical_lower_bound_complex_y;
					for (ii = 0; ii <= physical_upper_bound_complex_x; ii++)
					{
						other_image->complex_values[other_image->ReturnFourier1DAddressFromLogicalCoord(ii,jj,kk_logi)] = complex_values[ReturnFourier1DAddressFromLogicalCoord(ii,jj_logi,kk_logi)];
					}
				}


				// Deal with the positive Nyquist in the 3rd dimension
				kk = physical_index_of_first_negative_frequency_z;
				int kk_mirror = other_image->logical_z_dimension - physical_index_of_first_negative_frequency_z;
				//wxPrintf("\nkk = %i; kk_mirror = %i\n",kk,kk_mirror);
				int jj_mirror;
				//wxPrintf("Will loop jj from %i to %i\n",1,physical_index_of_first_negative_frequency_y);
				for (jj = 1; jj <= physical_index_of_first_negative_frequency_y; jj ++ )
				{
					//jj_mirror = other_image->logical_y_dimension - jj;
					jj_mirror = jj;
					for (ii = 0; ii <= physical_upper_bound_complex_x; ii++ )
					{
						//wxPrintf("(1) ii = %i; jj = %i; kk = %i; jj_mirror = %i; kk_mirror = %i\n",ii,jj,kk,jj_mirror,kk_mirror);
						other_image->complex_values[other_image-> ReturnFourier1DAddressFromPhysicalCoord(ii,jj,kk)] = other_image->complex_values[other_image->ReturnFourier1DAddressFromPhysicalCoord(ii,jj_mirror,kk_mirror)];
					}
				}
				//wxPrintf("Will loop jj from %i to %i\n", other_image->logical_y_dimension - physical_index_of_first_negative_frequency_y, other_image->logical_y_dimension - 1);
				for (jj = other_image->logical_y_dimension - physical_index_of_first_negative_frequency_y; jj <= other_image->logical_y_dimension - 1; jj ++)
				{
					//jj_mirror = other_image->logical_y_dimension - jj;
					jj_mirror = jj;
					for (ii = 0; ii <= physical_upper_bound_complex_x; ii++ )
					{
						//wxPrintf("(2) ii = %i; jj = %i; kk = %i; jj_mirror = %i; kk_mirror = %i\n",ii,jj,kk,jj_mirror,kk_mirror);
						other_image->complex_values[other_image-> ReturnFourier1DAddressFromPhysicalCoord(ii,jj,kk)] = other_image->complex_values[other_image->ReturnFourier1DAddressFromPhysicalCoord(ii,jj_mirror,kk_mirror)];
					}
				}
				jj = 0;
				for (ii = 0; ii <= physical_upper_bound_complex_x; ii++)
				{
					other_image->complex_values[other_image->ReturnFourier1DAddressFromPhysicalCoord(ii,jj,kk)] = other_image->complex_values[other_image->ReturnFourier1DAddressFromPhysicalCoord(ii,jj,kk_mirror)];
				}

			}
		}


	}

}


// Bilinear interpolation in real space, at point (x,y) where x and y are physical coordinates (i.e. first pixel has x,y = 0,0)
void Image::GetRealValueByLinearInterpolationNoBoundsCheckImage(float &x, float &y, float &interpolated_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(logical_z_dimension == 1, "Not for volumes");
	MyDebugAssertTrue(is_in_real_space, "Need to be in real space");

	const int i_start = int(x);
	const int j_start = int(y);
	const float x_dist = x - float(i_start);
	const float y_dist = y - float(j_start);
	const float x_dist_m = 1.0 - x_dist;
	const float y_dist_m = 1.0 - y_dist;

	const int address_1 = j_start * (logical_x_dimension + padding_jump_value) + i_start;
	const int address_2 = address_1 + logical_x_dimension + padding_jump_value;

	MyDebugAssertTrue(address_1+1 <= real_memory_allocated && address_1 >= 0,"Out of bounds, address 1\n");
	MyDebugAssertTrue(address_2+1 <= real_memory_allocated && address_2 >= 0,"Out of bounds, address 2\n");

	interpolated_value =    x_dist_m * y_dist_m * real_values[address_1]
						+	x_dist   * y_dist_m * real_values[address_1 + 1]
						+   x_dist_m * y_dist   * real_values[address_2]
						+   x_dist   * y_dist   * real_values[address_2 + 1];

}


void Image::Resize(int wanted_x_dimension, int wanted_y_dimension, int wanted_z_dimension, float wanted_padding_value)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(wanted_x_dimension != 0 && wanted_y_dimension != 0 && wanted_z_dimension != 0, "Resize dimension is zero");

	if (logical_x_dimension == wanted_x_dimension && logical_y_dimension == wanted_y_dimension && logical_z_dimension == wanted_z_dimension) return;

	Image temp_image;

	temp_image.Allocate(wanted_x_dimension, wanted_y_dimension, wanted_z_dimension, is_in_real_space);
	ClipInto(&temp_image, wanted_padding_value);

	//CopyFrom(&temp_image);
	Consume(&temp_image);
}


void Image::CopyFrom(Image *other_image)
{
	*this = other_image;
}

void Image::CopyLoopingAndAddressingFrom(Image *other_image)
{
	object_is_centred_in_box = other_image->object_is_centred_in_box;
	logical_x_dimension = other_image->logical_x_dimension;
	logical_y_dimension = other_image->logical_y_dimension;
	logical_z_dimension = other_image->logical_z_dimension;

	physical_upper_bound_complex_x = other_image->physical_upper_bound_complex_x;
	physical_upper_bound_complex_y = other_image->physical_upper_bound_complex_y;
	physical_upper_bound_complex_z = other_image->physical_upper_bound_complex_z;

	physical_address_of_box_center_x = other_image->physical_address_of_box_center_x;
	physical_address_of_box_center_y = other_image->physical_address_of_box_center_y;
	physical_address_of_box_center_z = other_image->physical_address_of_box_center_z;

	//physical_index_of_first_negative_frequency_x = other_image->physical_index_of_first_negative_frequency_x;
	physical_index_of_first_negative_frequency_y = other_image->physical_index_of_first_negative_frequency_y;
	physical_index_of_first_negative_frequency_z = other_image->physical_index_of_first_negative_frequency_z;

	fourier_voxel_size_x = other_image->fourier_voxel_size_x;
	fourier_voxel_size_y = other_image->fourier_voxel_size_y;
	fourier_voxel_size_z = other_image->fourier_voxel_size_z;

	logical_upper_bound_complex_x = other_image->logical_upper_bound_complex_x;
	logical_upper_bound_complex_y = other_image->logical_upper_bound_complex_y;
	logical_upper_bound_complex_z = other_image->logical_upper_bound_complex_z;

	logical_lower_bound_complex_x = other_image->logical_lower_bound_complex_x;
	logical_lower_bound_complex_y = other_image->logical_lower_bound_complex_y;
	logical_lower_bound_complex_z = other_image->logical_lower_bound_complex_z;

	logical_upper_bound_real_x = other_image->logical_upper_bound_complex_x;
	logical_upper_bound_real_y = other_image->logical_upper_bound_complex_y;
	logical_upper_bound_real_z = other_image->logical_upper_bound_complex_z;

	logical_lower_bound_real_x = other_image->logical_lower_bound_complex_x;
	logical_lower_bound_real_y = other_image->logical_lower_bound_complex_y;
	logical_lower_bound_real_z = other_image->logical_lower_bound_complex_z;

	padding_jump_value = other_image->padding_jump_value;
}

void Image::Consume(Image *other_image) // copy the parameters then directly steal the memory of another image, leaving it an empty shell
{
	MyDebugAssertTrue(other_image->is_in_memory, "Other image Memory not allocated");

	if (is_in_memory == true)
	{
		Deallocate();
	}

	is_in_real_space = other_image->is_in_real_space;
	real_memory_allocated = other_image->real_memory_allocated;
	CopyLoopingAndAddressingFrom(other_image);

	real_values = other_image->real_values;
	complex_values = other_image->complex_values;
	is_in_memory = other_image->is_in_memory;

	/* IMOD
  plan_fwd = other_image->plan_fwd;
	plan_bwd = other_image->plan_bwd;
	planned = other_image->planned;*/

	other_image->real_values = NULL;
	other_image->complex_values = NULL;
	other_image->is_in_memory = false;

	/* IMOD
  other_image->plan_fwd = NULL;
	other_image->plan_bwd = NULL;
	other_image->planned = false;*/

	number_of_real_space_pixels = other_image->number_of_real_space_pixels;
	ft_normalization_factor = other_image->ft_normalization_factor;

}

/* IMOD avoid remainderf
void Image::RealSpaceIntegerShift(int wanted_x_shift, int wanted_y_shift, int wanted_z_shift)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i, j, k;
	long pixel_counter = 0;
	long shifted_counter;
	float *buffer = new float[number_of_real_space_pixels];

	shifted_counter = - wanted_x_shift - logical_x_dimension * (wanted_y_shift + logical_y_dimension * wanted_z_shift);
	shifted_counter = remainderf(float(shifted_counter), float(number_of_real_space_pixels));
	if (shifted_counter < 0) shifted_counter += number_of_real_space_pixels;

	for (k = 0; k < logical_z_dimension; k++)
	{
		for (j = 0; j < logical_y_dimension; j++)
		{
			for (i = 0; i < logical_x_dimension; i++)
			{
				buffer[shifted_counter] = real_values[pixel_counter];
				pixel_counter++;
				shifted_counter++;
				if (shifted_counter >= number_of_real_space_pixels) shifted_counter -= number_of_real_space_pixels;
			}
			pixel_counter += padding_jump_value;
		}
	}
	shifted_counter = 0;
	pixel_counter = 0;
	for (k = 0; k < logical_z_dimension; k++)
	{
		for (j = 0; j < logical_y_dimension; j++)
		{
			for (i = 0; i < logical_x_dimension; i++)
			{
				real_values[pixel_counter] = buffer[shifted_counter];
				pixel_counter++;
				shifted_counter++;
			}
			pixel_counter += padding_jump_value;
		}
	}
	delete [] buffer;
}

void Image::DilateBinarizedMask(float dilation_radius)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Not in real space");

	int i, j, k;
	int l, m, n, m2, n2;
	int lim, nlim;
	long pixel_counter = 0;
	long shifted_counter = 0;
	float dilation_radius_squared = powf(dilation_radius, 2);
	float *buffer = new float[number_of_real_space_pixels];

	for (k = 0; k < logical_z_dimension; k++)
	{
		for (j = 0; j < logical_y_dimension; j++)
		{
			for (i = 0; i < logical_x_dimension; i++)
			{
				buffer[shifted_counter] = real_values[pixel_counter];
				pixel_counter++;
				shifted_counter++;
			}
			pixel_counter += padding_jump_value;
		}
	}

	lim = myroundint(dilation_radius);
	if (IsEven(lim)) lim++;
	nlim = lim;
	if (logical_z_dimension == 1) nlim = 0;
	for (n = -nlim; n <= nlim; n += 2)
	{
		n2 = n * n;
		for (m = -lim; m <= lim; m += 2)
		{
			m2 = m * m;
			for (l = -lim; l <= lim; l += 2)
			{
				if (l * l + m2 + n2 <= dilation_radius_squared)
				{
					shifted_counter = - l - logical_x_dimension * (m + logical_y_dimension * n);
					shifted_counter = remainderf(float(shifted_counter), float(number_of_real_space_pixels));
					if (shifted_counter < 0) shifted_counter += number_of_real_space_pixels;

					pixel_counter = 0;
					for (k = 0; k < logical_z_dimension; k++)
					{
						for (j = 0; j < logical_y_dimension; j++)
						{
							for (i = 0; i < logical_x_dimension; i++)
							{
								real_values[pixel_counter] += buffer[shifted_counter];
								pixel_counter++;
								shifted_counter++;
								if (shifted_counter >= number_of_real_space_pixels) shifted_counter -= number_of_real_space_pixels;
							}
							pixel_counter += padding_jump_value;
						}
					}
				}
			}
		}
	}

	pixel_counter = 0;
	for (k = 0; k < logical_z_dimension; k++)
	{
		for (j = 0; j < logical_y_dimension; j++)
		{
			for (i = 0; i < logical_x_dimension; i++)
			{
				if (real_values[pixel_counter] != 0.0) real_values[pixel_counter] = 1.0;
				pixel_counter++;
			}
			pixel_counter += padding_jump_value;
		}
	}
	delete [] buffer;
}*/

/* IMOD
void Image::PhaseShift(float wanted_x_shift, float wanted_y_shift, float wanted_z_shift)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");


	bool need_to_fft = false;

	long pixel_counter = 0;

	int k;
	int k_logical;

	int j;
	int j_logical;

	int i;

	float phase_z;
	float phase_y;
	float phase_x;

	std::complex<float> total_phase_shift;

	if (is_in_real_space == true)
	{
		ForwardFFT();
		need_to_fft = true;
	}

	for (k=0; k <= physical_upper_bound_complex_z; k++)
	{
		k_logical = ReturnFourierLogicalCoordGivenPhysicalCoord_Z(k);
		phase_z = ReturnPhaseFromShift(wanted_z_shift, k_logical, logical_z_dimension);

		for (j = 0; j <= physical_upper_bound_complex_y; j++)
		{
			j_logical = ReturnFourierLogicalCoordGivenPhysicalCoord_Y(j);
			phase_y = ReturnPhaseFromShift(wanted_y_shift, j_logical, logical_y_dimension);

			for (i = 0; i <= physical_upper_bound_complex_x; i++)
			{

				phase_x = ReturnPhaseFromShift(wanted_x_shift, i, logical_x_dimension);

				total_phase_shift = Return3DPhaseFromIndividualDimensions(phase_x, phase_y, phase_z);
				complex_values[pixel_counter] *= total_phase_shift;

				pixel_counter++;
			}
		}
	}

	if (need_to_fft == true) BackwardFFT();

}*/


bool Image::HasSameDimensionsAs(Image *other_image)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(other_image->is_in_memory, "Other image Memory not allocated");

	if (logical_x_dimension == other_image->logical_x_dimension && logical_y_dimension == other_image->logical_y_dimension && logical_z_dimension == other_image->logical_z_dimension) return true;
	else return false;
}

float Image::ReturnLinearInterpolated2D(float &wanted_physical_x_coordinate, float &wanted_physical_y_coordinate)
{
	MyDebugAssertTrue(is_in_memory, "Memory not allocated");
	MyDebugAssertTrue(is_in_real_space, "Is in Fourier space");

	if(wanted_physical_x_coordinate < 0 || wanted_physical_x_coordinate > logical_x_dimension - 1) return 0.0;
	if(wanted_physical_y_coordinate < 0 || wanted_physical_y_coordinate > logical_y_dimension - 1) return 0.0;

	int i;
	int j;
	int int_x_coordinate;
	int int_y_coordinate;
	int int_x_coordinate1;
	int int_y_coordinate1;
	int int_y;

	float weight_x;
	float weight_y;

	float sum = 0.0;

	int_x_coordinate = int(floorf(wanted_physical_x_coordinate));
	int_y_coordinate = int(floorf(wanted_physical_y_coordinate));
	int_x_coordinate1 = int_x_coordinate + 1;
	int_y_coordinate1 = int_y_coordinate + 1;
	int_x_coordinate1 = std::min(int_x_coordinate1, logical_x_dimension - 1);
	int_y_coordinate1 = std::min(int_y_coordinate1, logical_y_dimension - 1);

	for (j = int_y_coordinate; j <= int_y_coordinate1; j++)
	{
		weight_y = (1.0 - fabsf(wanted_physical_y_coordinate - j));
		int_y = (logical_x_dimension + padding_jump_value) * j;
		for (i = int_x_coordinate; i <= int_x_coordinate1; i++)
		{
			weight_x = (1.0 - fabsf(wanted_physical_x_coordinate - i));
			sum += real_values[int_y + i] * weight_x * weight_y;
		}
	}

	return sum;
}
void Image::CorrectMagnificationDistortion(float distortion_angle, float distortion_major_axis, float distortion_minor_axis)
{
	MyDebugAssertTrue(logical_z_dimension == 1, "Only 2D Images supported");

	long pixel_counter = 0;
	float angle_in_radians = deg_2_rad(distortion_angle);

	float x_scale_factor = 1.0 / distortion_major_axis;
	float y_scale_factor = 1.0 / distortion_minor_axis;

	float average_edge_value = ReturnAverageOfRealValuesOnEdges();

	float new_x;
	float new_y;

	float final_x;
	float final_y;

	int x,y;

	Image buffer_image;
	buffer_image.Allocate(logical_x_dimension, logical_y_dimension, logical_z_dimension, is_in_real_space);
	//buffer_image.CopyFrom(this);

	for (y = 0; y < logical_y_dimension; y++)
	{
		for (x = 0; x < logical_x_dimension; x++)
		{
			// first rotation

			new_x = float(y - physical_address_of_box_center_y) * sinf(-angle_in_radians) + float(x - physical_address_of_box_center_x) * cosf(-angle_in_radians);
			new_y = float(y - physical_address_of_box_center_y) * cosf(-angle_in_radians) - float(x - physical_address_of_box_center_x) * sinf(-angle_in_radians);

			// scale factor

			new_x *= x_scale_factor;
			new_y *= y_scale_factor;

			new_x += physical_address_of_box_center_x;
			new_y += physical_address_of_box_center_y;

			// rotate back

			final_x = float(new_y - physical_address_of_box_center_y) * sinf(angle_in_radians) + float(new_x - physical_address_of_box_center_x) * cosf(angle_in_radians);
			final_y = float(new_y - physical_address_of_box_center_y) * cosf(angle_in_radians) - float(new_x - physical_address_of_box_center_x) * sinf(angle_in_radians);

			final_x += physical_address_of_box_center_x;
			final_y += physical_address_of_box_center_y;

			if (final_x < 0 || final_x > logical_x_dimension - 1 || final_y < 0 || final_y > logical_y_dimension - 1) real_values[pixel_counter] = average_edge_value;
			else
			{
				buffer_image.real_values[pixel_counter] = ReturnLinearInterpolated2D(final_x, final_y);
			}

			pixel_counter++;
		}

		pixel_counter += padding_jump_value;
	}

	Consume(&buffer_image);
}
