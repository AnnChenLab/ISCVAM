import React, { useState } from 'react';
import { TextField, Button, Typography, Snackbar } from '@material-ui/core';

const Contact = () => {
    const [formData, setFormData] = useState({ name: '', email: '', message: '' });
    const [formErrors, setFormErrors] = useState({});
    const [successMessage, setSuccessMessage] = useState(false);

    const handleChange = (e) => {
        const { name, value } = e.target;
        setFormData({
            ...formData,
            [name]: value,
        });
    };

    const validateForm = () => {
        const errors = {};
        if (!formData.name) errors.name = 'Name is required';
        if (!formData.email) errors.email = 'Email is required';
        else if (!/\S+@\S+\.\S+/.test(formData.email)) errors.email = 'Email address is invalid';
        if (!formData.message) errors.message = 'Message is required';
        return errors;
    };

    const handleSubmit = (e) => {
        e.preventDefault();
        const errors = validateForm();
        if (Object.keys(errors).length === 0) {
            setFormErrors({});
            setSuccessMessage(true);
            // Clear the form after successful submission
            setFormData({ name: '', email: '', message: '' });
            // You can add form submission logic here (e.g., API calls)
        } else {
            setFormErrors(errors);
        }
    };

    return (
        <div style={{ maxWidth: '600px', margin: '0 auto', padding: '20px' }}>
            <form onSubmit={handleSubmit} noValidate>
                <Typography variant="h5" gutterBottom>
                    Contact Us
                </Typography>

                <TextField
                    label="Name"
                    name="name"
                    fullWidth
                    variant="outlined"
                    value={formData.name}
                    onChange={handleChange}
                    error={!!formErrors.name}
                    helperText={formErrors.name}
                    margin="normal"
                />

                <TextField
                    label="Email"
                    name="email"
                    type="email"
                    fullWidth
                    variant="outlined"
                    value={formData.email}
                    onChange={handleChange}
                    error={!!formErrors.email}
                    helperText={formErrors.email}
                    margin="normal"
                />

                <TextField
                    label="Message"
                    name="message"
                    fullWidth
                    variant="outlined"
                    multiline
                    rows={4}
                    value={formData.message}
                    onChange={handleChange}
                    error={!!formErrors.message}
                    helperText={formErrors.message}
                    margin="normal"
                />

                <Button type="submit" variant="contained" color="primary" style={{ marginTop: '20px' }}>
                    Send Message
                </Button>

                {/* Success Snackbar */}
                <Snackbar
                    open={successMessage}
                    autoHideDuration={4000}
                    onClose={() => setSuccessMessage(false)}
                    message="Message sent successfully!"
                />
            </form>
        </div>
    );
};

export default Contact;
