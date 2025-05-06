import React, { useState, useEffect } from 'react';
import { TextField, Button, Typography, Snackbar } from '@material-ui/core';
import axios from 'axios';

const Contact = () => {
    const [formData, setFormData] = useState({ name: '', email: '', message: '' });
    const [formErrors, setFormErrors] = useState({});
    const [successMessage, setSuccessMessage] = useState(false);
    const [errorMessage, setErrorMessage] = useState(false);
    const [messages, setMessages] = useState([]);

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

    const handleSubmit = async (e) => {
        e.preventDefault();
        const errors = validateForm();
        if (Object.keys(errors).length === 0) {
            setFormErrors({});
            try {
                await axios.post('https://chenlab.chpc.utah.edu/iscvam/backend/api/save-message', formData);
                setSuccessMessage(true);
                setFormData({ name: '', email: '', message: '' });
                fetchMessages(); // Fetch updated messages after a new one is saved
            } catch (error) {
                console.error('Error saving message:', error);
                setErrorMessage(true);
            }
        } else {
            setFormErrors(errors);
        }
    };

    const fetchMessages = async () => {
        try {
            const response = await axios.get('https://chenlab.chpc.utah.edu/iscvam/backend/api/messages');
            setMessages(response.data);
        } catch (error) {
            console.error('Error fetching messages:', error);
        }
    };

    useEffect(() => {
        fetchMessages();
    }, []);

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
                    message="Message saved successfully!"
                />

                {/* Error Snackbar */}
                <Snackbar
                    open={errorMessage}
                    autoHideDuration={4000}
                    onClose={() => setErrorMessage(false)}
                    message="Failed to save the message."
                />
            </form>
        </div>
    );
};

export default Contact;
