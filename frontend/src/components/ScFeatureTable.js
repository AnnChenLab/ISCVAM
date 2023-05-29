import React, { useState } from "react";
import { useTheme, Table, TableBody, TableCell, TableContainer, 
         TableFooter, TablePagination, TableRow, Box, IconButton, 
         Icon } from '@mui/material';
import PropTypes from 'prop-types';

function TablePaginationActions({ count, page, rowsPerPage, onPageChange }) {
    const theme = useTheme();
    const isRtl = theme.direction === 'rtl';

    const handleButtonClick = (newPage) => (event) => onPageChange(event, newPage);

    return (
        <Box sx={{ flexShrink: 0, ml: 2.5 }}>
            <IconButton onClick={handleButtonClick(0)} disabled={page === 0} aria-label="first page">
                <Icon>{isRtl ? 'last_page' : 'first_page'}</Icon>
            </IconButton>
            <IconButton onClick={handleButtonClick(page - 1)} disabled={page === 0} aria-label="previous page">
                <Icon>{isRtl ? 'keyboard_arrow_right' : 'keyboard_arrow_left'}</Icon>
            </IconButton>
            <IconButton onClick={handleButtonClick(page + 1)} disabled={page >= Math.ceil(count / rowsPerPage) - 1} aria-label="next page">
                <Icon>{isRtl ? 'keyboard_arrow_left' : 'keyboard_arrow_right'}</Icon>
            </IconButton>
            <IconButton onClick={handleButtonClick(Math.max(0, Math.ceil(count / rowsPerPage) - 1))} disabled={page >= Math.ceil(count / rowsPerPage) - 1} aria-label="last page">
                <Icon>{isRtl ? 'first_page' : 'last_page'}</Icon>
            </IconButton>
        </Box>
    );
}

TablePaginationActions.propTypes = {
    count: PropTypes.number.isRequired,
    onPageChange: PropTypes.func.isRequired,
    page: PropTypes.number.isRequired,
    rowsPerPage: PropTypes.number.isRequired,
};

export function FeatureTable({features, assay}) {
    const [page, setPage] = useState(0);
    const [rowsPerPage, setRowsPerPage] = useState(50);
    const colsPerRow = 2;

    const handleChangePage = (_, newPage) => setPage(newPage);

    const handleChangeRowsPerPage = (event) => {
        setRowsPerPage(parseInt(event.target.value, 10));
        setPage(0);
    };

    return (
        <TableContainer>
            <Table sx={{ maxWidth: 600, marginTop: 2, captionSide: 'top' }} size="small" aria-label="features table">
                {/* <caption>Available {assay} features</caption> */}
                <TableBody>
                    {[...Array(rowsPerPage).keys()].map((rowid) => (
                        <TableRow key={rowid}>
                            <TableCell component="th" scope="row">
                                {features[page * rowsPerPage * colsPerRow + rowid * colsPerRow]}
                            </TableCell>
                            {[...Array(colsPerRow - 1).keys()].map(cidx => (
                                <TableCell align="right" key={cidx}>
                                    {features[page * rowsPerPage * colsPerRow + rowid * colsPerRow + 1 + cidx] || ''}
                                </TableCell>
                            ))}
                        </TableRow>
                    ))}
                </TableBody>
                <TableFooter>
                    <TableRow>
                        <TablePagination
                            rowsPerPageOptions={[5, 15, 50]}
                            colSpan={3}
                            count={Math.ceil(features.length / colsPerRow)}
                            rowsPerPage={rowsPerPage}
                            page={page}
                            SelectProps={{
                                inputProps: { 'aria-label': 'rows per page' },
                                native: true,
                            }}
                            onPageChange={handleChangePage}
                            onRowsPerPageChange={handleChangeRowsPerPage}
                            ActionsComponent={TablePaginationActions}
                        />
                    </TableRow>
                </TableFooter>
            </Table>
        </TableContainer>
    );
}
